#define kiss_fft_scalar double
#include "stft_fft_plugin.h"
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include "signal_helpers.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 1D 리샘플 — 다운샘플(n_src > n_dst)은 max-pooling(peak 보존), 업샘플/동일은
 * 선형 보간. 오더 라인 같은 가느다란 능선이 다운샘플 시 평균화로 사라지는 것을
 * 막기 위해, 축소 방향에서는 source span 의 최댓값을 취한다.
 * (입력은 dB 도메인 값이므로 max 가 곧 peak 보존) */
static void resample_1d(const double* src, int n_src, double* dst, int n_dst) {
    if (n_dst <= 0) return;
    if (n_src <= 0) {
        for (int i = 0; i < n_dst; i++) dst[i] = 0.0;
        return;
    }
    if (n_src == 1) {
        for (int i = 0; i < n_dst; i++) dst[i] = src[0];
        return;
    }
    if (n_dst == 1) {
        /* 단일 출력: 전체 구간 최댓값(peak 보존) */
        double m = src[0];
        for (int i = 1; i < n_src; i++) if (src[i] > m) m = src[i];
        dst[0] = m;
        return;
    }
    if (n_dst < n_src) {
        /* 다운샘플: 각 출력 cell 이 덮는 source span 의 max */
        for (int i = 0; i < n_dst; i++) {
            int s0 = (int)((long long)i * n_src / n_dst);
            int s1 = (int)((long long)(i + 1) * n_src / n_dst);
            if (s1 <= s0) s1 = s0 + 1;
            if (s1 > n_src) s1 = n_src;
            double m = src[s0];
            for (int s = s0 + 1; s < s1; s++) if (src[s] > m) m = src[s];
            dst[i] = m;
        }
    } else {
        /* 업샘플/동일: 선형 보간 (n_dst==n_src 이면 ratio=1 로 그대로 복사) */
        double ratio = (double)(n_src - 1) / (double)(n_dst - 1);
        for (int i = 0; i < n_dst; i++) {
            double sp = i * ratio;
            int s0 = (int)sp;
            int s1 = s0 + 1;
            if (s1 >= n_src) s1 = n_src - 1;
            double w = sp - s0;
            dst[i] = src[s0] * (1.0 - w) + src[s1] * w;
        }
    }
}

/* 원본 STFT grid(num_frames × num_bins)에서 [frame_lo,frame_hi) × [bin_lo,bin_hi)
 * 대역만 잘라(crop) target_width × target_height 로 리샘플한다(crop-then-resample).
 * 주파수(열)·시간(행) 축을 분리해 각각 resample_1d 적용 — 축소 축은 max-pooling.
 * dst 레이아웃: dst[w * target_height + h] (= 기존 resampled 와 동일).
 * 성공 0, 실패 -1. */
static int resample_crop_region(
    const double* spec, int num_bins,
    int frame_lo, int frame_hi, int bin_lo, int bin_hi,
    double* dst, int target_width, int target_height
) {
    int crop_rows = frame_hi - frame_lo;   /* 시간 프레임 수 */
    int crop_cols = bin_hi - bin_lo;       /* 주파수 bin 수 */
    if (crop_rows <= 0 || crop_cols <= 0 || target_width <= 0 || target_height <= 0) {
        return -1;
    }

    /* 중간 버퍼: 주파수축만 먼저 리샘플 → tmp[crop_rows][target_height] */
    double* tmp = (double*)malloc((size_t)crop_rows * target_height * sizeof(double));
    double* col_src = (double*)malloc((size_t)crop_cols * sizeof(double));
    double* col_dst = (double*)malloc((size_t)target_height * sizeof(double));
    if (!tmp || !col_src || !col_dst) { free(tmp); free(col_src); free(col_dst); return -1; }

    for (int r = 0; r < crop_rows; r++) {
        const double* row = spec + (size_t)(frame_lo + r) * num_bins + bin_lo;
        for (int c = 0; c < crop_cols; c++) col_src[c] = row[c];
        resample_1d(col_src, crop_cols, col_dst, target_height);
        for (int h = 0; h < target_height; h++) {
            tmp[(size_t)r * target_height + h] = col_dst[h];
        }
    }
    free(col_src);

    /* 시간축 리샘플: 각 출력 주파수열에 대해 crop_rows → target_width */
    double* row_src = (double*)malloc((size_t)crop_rows * sizeof(double));
    double* row_dst = (double*)malloc((size_t)target_width * sizeof(double));
    if (!row_src || !row_dst) { free(tmp); free(col_dst); free(row_src); free(row_dst); return -1; }

    for (int h = 0; h < target_height; h++) {
        for (int r = 0; r < crop_rows; r++) row_src[r] = tmp[(size_t)r * target_height + h];
        resample_1d(row_src, crop_rows, row_dst, target_width);
        for (int w = 0; w < target_width; w++) {
            dst[(size_t)w * target_height + h] = row_dst[w];
        }
    }

    free(tmp); free(col_dst); free(row_src); free(row_dst);
    return 0;
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
    int32_t weighting_type,
    double overlap_percent,
    double rms_range_multiplier,
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

    /* 3. STFT: frame-by-frame FFT.
     * hop_size 는 overlap_percent 로부터 계산. 기본값 50 일 때 chunk_size/2
     * 와 byte-exact 동등 (chunk_size 가 nextPow2 이므로 짝수). */
    int hop_size = (int)((double)chunk_size * (1.0 - overlap_percent / 100.0) + 0.5);
    if (hop_size < 1) hop_size = 1;
    if (hop_size > chunk_size) hop_size = chunk_size;
    int num_bins = chunk_size / 2 + 1;
    int num_frames = 0;
    for (int s = 0; s + chunk_size <= padded_len; s += hop_size) num_frames++;
    if (num_frames == 0) { free(padded); free(window); return NULL; }

    double* windowed = (double*)malloc(chunk_size * sizeof(double));
    kiss_fft_cpx* fft_out = (kiss_fft_cpx*)malloc(num_bins * sizeof(kiss_fft_cpx));
    double* spectrogram = (double*)malloc(num_frames * num_bins * sizeof(double));
    /* [이슈 #2 / ADR-003] dB 변환 직전 단계의 선형 amplitude 보관. */
    double* amplitude_raw = (double*)malloc(num_frames * num_bins * sizeof(double));

    if (!windowed || !fft_out || !spectrogram || !amplitude_raw) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram); free(amplitude_raw);
        return NULL;
    }

    kiss_fftr_cfg cfg = kiss_fftr_alloc(chunk_size, 0, NULL, NULL);
    if (!cfg) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram); free(amplitude_raw);
        return NULL;
    }

    double inv_chunk = 1.0 / (double)chunk_size;
    double ln10 = log(10.0);
    double factor = 20.0 / ln10;
    double db_offset_val = is_acc
        ? factor * log(9.807 / db_ref)
        : factor * log(1.0 / db_ref);
    double sqrt2 = sqrt(2.0);

    /* amplitude 임시 버퍼 — weighting/smoothing 을 적용한 뒤 dB 변환 */
    double* frame_amp = (double*)malloc(num_bins * sizeof(double));
    if (!frame_amp) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram); free(amplitude_raw);
        kiss_fft_free(cfg);
        return NULL;
    }

    /* 주파수 가중치 prefactor: 10^(W(f)/20). linear (type=0) 일 때 1.0 (no-op). */
    double* weight_factor = (double*)malloc(num_bins * sizeof(double));
    if (!weight_factor) {
        free(padded); free(window); free(windowed); free(fft_out); free(spectrogram); free(amplitude_raw); free(frame_amp);
        kiss_fft_free(cfg);
        return NULL;
    }
    if (weighting_type > 0) {
        double df = sampling_frequency / (double)chunk_size;
        for (int i = 0; i < num_bins; i++) {
            double freq = (double)i * df;
            double w_db = get_weighting_db(freq, weighting_type);
            weight_factor[i] = pow(10.0, w_db / 20.0);
        }
    } else {
        for (int i = 0; i < num_bins; i++) weight_factor[i] = 1.0;
    }

    int smooth_bins = (int)rms_range_multiplier;

    int frame_idx = 0;
    for (int s = 0; s + chunk_size <= padded_len; s += hop_size) {
        for (int i = 0; i < chunk_size; i++) {
            windowed[i] = padded[s + i] * window[i];
        }
        kiss_fftr(cfg, windowed, fft_out);

        /* 3a. amplitude 계산 (RMS-scaled) + weighting */
        for (int i = 0; i < num_bins; i++) {
            double amp = sqrt(fft_out[i].r * fft_out[i].r + fft_out[i].i * fft_out[i].i);
            double norm = (i != 0 && i != num_bins - 1)
                ? amp * inv_chunk * 2.0
                : amp * inv_chunk;
            frame_amp[i] = norm * sqrt2 * weight_factor[i];
        }

        /* 3b. RMS smoothing (smooth_bins<=1 이면 no-op) */
        apply_rms_smoothing(frame_amp, num_bins, smooth_bins);

        /* 3c. dB 변환 + amplitude_raw 동시 저장 (이슈 #2: dB 변환 직전 선형 값) */
        for (int i = 0; i < num_bins; i++) {
            double x = frame_amp[i];
            double db = (x > 0.0) ? db_offset_val + factor * log(x) : -200.0;
            spectrogram[frame_idx * num_bins + i] = db;
            amplitude_raw[frame_idx * num_bins + i] = x;
        }
        frame_idx++;
    }

    free(padded);
    free(window);
    free(windowed);
    free(fft_out);
    free(frame_amp);
    free(weight_factor);
    kiss_fft_free(cfg);

    /* 4. Crop-then-resample (STFT 줌 시 주파수 해상도 미회복 수정).
     *
     * 기존엔 원본 grid 를 곧바로 target 크기로 다운샘플한 뒤 crop 했다
     * (resample→crop). 그 결과 화면 주파수 해상도가 Nyquist/target_height 로
     * 고정되어, 대역을 좁혀도(줌) 해상도가 회복되지 않았다.
     *
     * 이제 원본 grid(num_frames × num_bins)에서 선택 대역을 먼저 crop 한 뒤,
     * 그 대역만 target_width × target_height 로 리샘플한다(crop→resample).
     * 다운 비율>1 축은 max-pooling 으로 peak(오더 능선) 보존.
     *
     * [이슈 #2] amplitude 는 여전히 raw grid(num_frames × num_bins) 그대로 노출. */
    double df = sampling_frequency / (double)chunk_size;   /* raw bin 폭 */
    if (df <= 0.0) df = 1.0;

    /* 주파수 → raw bin 범위 [bin_lo, bin_hi) */
    int bin_lo = (int)(start_freq_hz / df);
    int bin_hi = (int)ceil(end_freq_hz / df);
    if (bin_lo < 0) bin_lo = 0;
    if (bin_hi > num_bins) bin_hi = num_bins;
    if (bin_hi <= bin_lo) bin_hi = (bin_lo < num_bins) ? bin_lo + 1 : num_bins;
    if (bin_lo >= bin_hi) bin_lo = bin_hi - 1;

    /* 시간 → raw frame 범위 [frame_lo, frame_hi) */
    int frame_lo = (total_duration > 0) ? (int)(start_time / total_duration * num_frames) : 0;
    int frame_hi = (total_duration > 0) ? (int)ceil(end_time / total_duration * num_frames) : num_frames;
    if (frame_lo < 0) frame_lo = 0;
    if (frame_hi > num_frames) frame_hi = num_frames;
    if (frame_hi <= frame_lo) frame_hi = (frame_lo < num_frames) ? frame_lo + 1 : num_frames;
    if (frame_lo >= frame_hi) frame_lo = frame_hi - 1;

    /* 선택 대역을 target 크기로 리샘플. 출력이 곧 화면 픽셀 grid 와 동일하므로
     * 별도 crop 단계가 필요 없다 (대역이 전체 target_height 를 채움 = 해상도 회복). */
    double* resampled = (double*)malloc((size_t)target_width * target_height * sizeof(double));
    if (!resampled) {
        free(spectrogram); free(amplitude_raw);
        return NULL;
    }
    if (resample_crop_region(spectrogram, num_bins,
                             frame_lo, frame_hi, bin_lo, bin_hi,
                             resampled, target_width, target_height) != 0) {
        free(spectrogram); free(amplitude_raw); free(resampled);
        return NULL;
    }

    int filtered_width = target_width;
    int filtered_height = target_height;

    /* 5. Auto min/max — crop 된 고해상도 대역(resampled) 위에서 계산 */
    double min_val = min_level, max_val = max_level;
    if (auto_calc_min_max) {
        min_val = 1e30;
        max_val = -1e30;
        for (int i = 0; i < filtered_width * filtered_height; i++) {
            double v = resampled[i];
            if (isfinite(v)) {
                if (v < min_val) min_val = v;
                if (v > max_val) max_val = v;
            }
        }
        if (min_val > max_val) { min_val = 0; max_val = 1; }
    }

    /* 6. Generate pixels (resampled 가 곧 표시 대역) */
    uint8_t* pixels = (uint8_t*)malloc((size_t)filtered_width * filtered_height * 4);
    if (!pixels) { free(spectrogram); free(amplitude_raw); free(resampled); return NULL; }

    generate_pixels_from_data(resampled, filtered_width, filtered_height,
                              color_lut, color_lut_size, min_val, max_val, pixels);

    /* Build result */
    StftResult* result = (StftResult*)malloc(sizeof(StftResult));
    if (!result) { free(spectrogram); free(amplitude_raw); free(resampled); free(pixels); return NULL; }

    result->spectrogram = resampled;
    result->spec_rows = target_width;
    result->spec_cols = target_height;
    /* [이슈 #2 수정] amplitude 는 raw grid (num_frames × num_bins) 로 노출. */
    result->amplitude = amplitude_raw;
    result->amp_rows = num_frames;
    result->amp_cols = num_bins;
    result->pixels = pixels;
    result->pixel_width = filtered_width;
    result->pixel_height = filtered_height;
    result->min_level = min_val;
    result->max_level = max_val;

    free(spectrogram);
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
    free(result->amplitude);
    free(result->pixels);
    free(result);
}
