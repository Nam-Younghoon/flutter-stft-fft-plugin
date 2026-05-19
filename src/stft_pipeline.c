#define kiss_fft_scalar double
#include "stft_fft_plugin.h"
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void bilinear_resample(
    const double* src, int src_rows, int src_cols,
    double* dst, int dst_rows, int dst_cols
) {
    if (dst_rows <= 1 || dst_cols <= 1 || src_rows <= 1 || src_cols <= 1) {
        memset(dst, 0, dst_rows * dst_cols * sizeof(double));
        return;
    }

    double row_ratio = (double)(src_rows - 1) / (double)(dst_rows - 1);
    double col_ratio = (double)(src_cols - 1) / (double)(dst_cols - 1);

    for (int r = 0; r < dst_rows; r++) {
        double src_r = r * row_ratio;
        int r0 = (int)src_r;
        int r1 = r0 + 1;
        if (r1 >= src_rows) r1 = src_rows - 1;
        double rw = src_r - r0;

        for (int c = 0; c < dst_cols; c++) {
            double src_c = c * col_ratio;
            int c0 = (int)src_c;
            int c1 = c0 + 1;
            if (c1 >= src_cols) c1 = src_cols - 1;
            double cw = src_c - c0;

            double v00 = src[r0 * src_cols + c0];
            double v01 = src[r0 * src_cols + c1];
            double v10 = src[r1 * src_cols + c0];
            double v11 = src[r1 * src_cols + c1];

            double v0 = v00 * (1.0 - cw) + v01 * cw;
            double v1 = v10 * (1.0 - cw) + v11 * cw;
            dst[r * dst_cols + c] = v0 * (1.0 - rw) + v1 * rw;
        }
    }
}

static void generate_pixels_from_data(
    const double* data, int width, int height,
    const int32_t* color_lut, int32_t color_lut_size,
    double min_level, double max_level,
    uint8_t* pixels
) {
    double level_range = max_level - min_level;
    if (level_range == 0.0) level_range = 1.0;
    double inv_range = 255.0 / level_range;
    int lut_max = color_lut_size - 1;

    for (int w = 0; w < width; w++) {
        for (int h = 0; h < height; h++) {
            double val = data[w * height + h];
            double normalized = (val - min_level) * inv_range;
            int idx = (int)normalized;
            if (idx < 0) idx = 0;
            if (idx > lut_max) idx = lut_max;

            int32_t packed = color_lut[idx];
            int pixel_index = (width * (height - 1 - h) + w) * 4;

            pixels[pixel_index + 0] = (packed >> 24) & 0xFF;
            pixels[pixel_index + 1] = (packed >> 16) & 0xFF;
            pixels[pixel_index + 2] = (packed >> 8) & 0xFF;
            pixels[pixel_index + 3] = packed & 0xFF;
        }
    }
}

FFI_PLUGIN_EXPORT StftResult* compute_stft_spectrogram(
    const double* input,
    int32_t input_length,
    int32_t chunk_size,
    double sampling_frequency,
    int32_t is_acc,
    double db_ref,
    int32_t target_width,
    int32_t target_height,
    double start_freq_hz,
    double end_freq_hz,
    double start_time,
    double end_time,
    double total_duration,
    const int32_t* color_lut,
    int32_t color_lut_size,
    int32_t auto_calc_min_max,
    double min_level,
    double max_level
) {
    if (!input || input_length < chunk_size || chunk_size < 2) return NULL;

    /* 1. Zero padding */
    int offset = chunk_size / 2;
    int padded_len = offset + input_length;
    double* padded = (double*)calloc(padded_len, sizeof(double));
    if (!padded) return NULL;
    memcpy(padded + offset, input, input_length * sizeof(double));

    /* 2. Hanning window */
    double* window = (double*)malloc(chunk_size * sizeof(double));
    if (!window) { free(padded); return NULL; }
    generate_hanning_window(window, chunk_size);

    /* 3. STFT: frame-by-frame FFT */
    int hop_size = chunk_size / 2;
    int num_bins = chunk_size / 2 + 1;
    int num_frames = 0;
    for (int s = 0; s + chunk_size <= padded_len; s += hop_size) num_frames++;
    if (num_frames == 0) { free(padded); free(window); return NULL; }

    double* windowed = (double*)malloc(chunk_size * sizeof(double));
    kiss_fft_cpx* fft_out = (kiss_fft_cpx*)malloc(num_bins * sizeof(kiss_fft_cpx));
    double* spectrogram = (double*)malloc(num_frames * num_bins * sizeof(double));

    if (!windowed || !fft_out || !spectrogram) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram);
        return NULL;
    }

    kiss_fftr_cfg cfg = kiss_fftr_alloc(chunk_size, 0, NULL, NULL);
    if (!cfg) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram);
        return NULL;
    }

    double inv_chunk = 1.0 / (double)chunk_size;
    double ln10 = log(10.0);
    double factor = 20.0 / ln10;
    double db_offset_val = is_acc
        ? factor * log(9.807 / db_ref)
        : factor * log(1.0 / db_ref);
    double sqrt2 = sqrt(2.0);

    int frame_idx = 0;
    for (int s = 0; s + chunk_size <= padded_len; s += hop_size) {
        for (int i = 0; i < chunk_size; i++) {
            windowed[i] = padded[s + i] * window[i];
        }
        kiss_fftr(cfg, windowed, fft_out);

        for (int i = 0; i < num_bins; i++) {
            double amp = sqrt(fft_out[i].r * fft_out[i].r + fft_out[i].i * fft_out[i].i);
            double norm = (i != 0 && i != num_bins - 1)
                ? amp * inv_chunk * 2.0
                : amp * inv_chunk;
            double x = norm * sqrt2;
            double db = (x > 0.0) ? db_offset_val + factor * log(x) : -200.0;
            spectrogram[frame_idx * num_bins + i] = db;
        }
        frame_idx++;
    }

    free(padded);
    free(window);
    free(windowed);
    free(fft_out);
    kiss_fft_free(cfg);

    /* 4. Bilinear resample to (target_width, target_height) */
    double* resampled = (double*)malloc(target_width * target_height * sizeof(double));
    if (!resampled) { free(spectrogram); return NULL; }
    bilinear_resample(spectrogram, num_frames, num_bins, resampled, target_width, target_height);

    /* 5. Frequency filtering */
    double nyquist = sampling_frequency / 2.0;
    int start_bin = (int)((start_freq_hz * target_height) / nyquist);
    int end_bin = (int)ceil((end_freq_hz * target_height) / nyquist);
    if (start_bin < 0) start_bin = 0;
    if (end_bin > target_height) end_bin = target_height;
    if (end_bin <= start_bin && start_bin < target_height) end_bin = start_bin + 1;
    int filtered_height = end_bin - start_bin;

    /* 6. Time filtering */
    int start_frame = (total_duration > 0) ? (int)(start_time / total_duration * target_width) : 0;
    int end_frame = (total_duration > 0) ? (int)ceil(end_time / total_duration * target_width) : target_width;
    if (start_frame < 0) start_frame = 0;
    if (end_frame > target_width) end_frame = target_width;
    if (end_frame <= start_frame) end_frame = start_frame + 1;
    int filtered_width = end_frame - start_frame;

    /* Extract filtered region */
    double* filtered = (double*)malloc(filtered_width * filtered_height * sizeof(double));
    if (!filtered) { free(spectrogram); free(resampled); return NULL; }

    for (int w = 0; w < filtered_width; w++) {
        for (int h = 0; h < filtered_height; h++) {
            filtered[w * filtered_height + h] = resampled[(start_frame + w) * target_height + (start_bin + h)];
        }
    }

    /* 7. Auto min/max */
    double min_val = min_level, max_val = max_level;
    if (auto_calc_min_max) {
        min_val = 1e30;
        max_val = -1e30;
        for (int i = 0; i < filtered_width * filtered_height; i++) {
            double v = filtered[i];
            if (isfinite(v)) {
                if (v < min_val) min_val = v;
                if (v > max_val) max_val = v;
            }
        }
        if (min_val > max_val) { min_val = 0; max_val = 1; }
    }

    /* 8. Generate pixels */
    uint8_t* pixels = (uint8_t*)malloc(filtered_width * filtered_height * 4);
    if (!pixels) { free(spectrogram); free(resampled); free(filtered); return NULL; }

    generate_pixels_from_data(filtered, filtered_width, filtered_height,
                              color_lut, color_lut_size, min_val, max_val, pixels);

    /* Build result */
    StftResult* result = (StftResult*)malloc(sizeof(StftResult));
    if (!result) { free(spectrogram); free(resampled); free(filtered); free(pixels); return NULL; }

    result->spectrogram = resampled;
    result->spec_rows = target_width;
    result->spec_cols = target_height;
    result->pixels = pixels;
    result->pixel_width = filtered_width;
    result->pixel_height = filtered_height;
    result->min_level = min_val;
    result->max_level = max_val;

    free(spectrogram);
    free(filtered);
    return result;
}

FFI_PLUGIN_EXPORT FilteredPixelResult* filter_and_generate_pixels(
    const double* spectrogram,
    int32_t spec_rows,
    int32_t spec_cols,
    double sampling_frequency,
    int32_t chunk_size,
    double total_duration,
    double start_freq_hz,
    double end_freq_hz,
    double start_time,
    double end_time,
    int32_t target_width,
    int32_t target_height,
    const int32_t* color_lut,
    int32_t color_lut_size,
    double min_level,
    double max_level
) {
    if (!spectrogram || spec_rows <= 0 || spec_cols <= 0) return NULL;
    (void)chunk_size; (void)target_width; (void)target_height;

    double nyquist = sampling_frequency / 2.0;
    int start_bin = (int)((start_freq_hz * spec_cols) / nyquist);
    int end_bin = (int)ceil((end_freq_hz * spec_cols) / nyquist);
    if (start_bin < 0) start_bin = 0;
    if (end_bin > spec_cols) end_bin = spec_cols;
    if (end_bin <= start_bin && start_bin < spec_cols) end_bin = start_bin + 1;
    int filtered_height = end_bin - start_bin;

    int start_frame = (total_duration > 0) ? (int)(start_time / total_duration * spec_rows) : 0;
    int end_frame = (total_duration > 0) ? (int)ceil(end_time / total_duration * spec_rows) : spec_rows;
    if (start_frame < 0) start_frame = 0;
    if (end_frame > spec_rows) end_frame = spec_rows;
    if (end_frame <= start_frame) end_frame = start_frame + 1;
    int filtered_width = end_frame - start_frame;

    double* filtered = (double*)malloc(filtered_width * filtered_height * sizeof(double));
    if (!filtered) return NULL;

    for (int w = 0; w < filtered_width; w++) {
        for (int h = 0; h < filtered_height; h++) {
            filtered[w * filtered_height + h] = spectrogram[(start_frame + w) * spec_cols + (start_bin + h)];
        }
    }

    uint8_t* pixels = (uint8_t*)malloc(filtered_width * filtered_height * 4);
    if (!pixels) { free(filtered); return NULL; }

    generate_pixels_from_data(filtered, filtered_width, filtered_height,
                              color_lut, color_lut_size, min_level, max_level, pixels);

    FilteredPixelResult* result = (FilteredPixelResult*)malloc(sizeof(FilteredPixelResult));
    if (!result) { free(filtered); free(pixels); return NULL; }

    result->filtered_data = filtered;
    result->filtered_rows = filtered_width;
    result->filtered_cols = filtered_height;
    result->pixels = pixels;
    result->pixel_width = filtered_width;
    result->pixel_height = filtered_height;
    return result;
}

FFI_PLUGIN_EXPORT void free_filtered_pixel_result(FilteredPixelResult* result) {
    if (!result) return;
    free(result->filtered_data);
    free(result->pixels);
    free(result);
}

FFI_PLUGIN_EXPORT PixelBufferResult* generate_pixel_buffer(
    const double* filtered_data,
    int32_t filtered_rows,
    int32_t filtered_cols,
    int32_t target_width,
    int32_t target_height,
    const int32_t* color_lut,
    int32_t color_lut_size,
    double min_level,
    double max_level
) {
    if (!filtered_data || filtered_rows <= 0 || filtered_cols <= 0) return NULL;
    (void)target_width; (void)target_height;

    int width = filtered_rows;
    int height = filtered_cols;

    uint8_t* pixels = (uint8_t*)malloc(width * height * 4);
    if (!pixels) return NULL;

    generate_pixels_from_data(filtered_data, width, height,
                              color_lut, color_lut_size, min_level, max_level, pixels);

    PixelBufferResult* result = (PixelBufferResult*)malloc(sizeof(PixelBufferResult));
    if (!result) { free(pixels); return NULL; }

    result->pixels = pixels;
    result->pixel_width = width;
    result->pixel_height = height;
    return result;
}

FFI_PLUGIN_EXPORT void free_pixel_buffer_result(PixelBufferResult* result) {
    if (!result) return;
    free(result->pixels);
    free(result);
}

FFI_PLUGIN_EXPORT void free_stft_result(StftResult* result) {
    if (!result) return;
    free(result->spectrogram);
    free(result->pixels);
    free(result);
}
