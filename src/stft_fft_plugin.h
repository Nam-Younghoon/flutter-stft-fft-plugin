#ifndef STFT_FFT_PLUGIN_H
#define STFT_FFT_PLUGIN_H

#include <stdint.h>
#include <stdlib.h>

#if _WIN32
#define FFI_PLUGIN_EXPORT __declspec(dllexport)
#else
#define FFI_PLUGIN_EXPORT __attribute__((visibility("default"))) __attribute__((used))
#endif

#ifdef __cplusplus
extern "C" {
#endif

FFI_PLUGIN_EXPORT int stft_fft_plugin_version(void);

/* ── FFT Engine (fft_engine.c) ── */

typedef struct {
    double* frequencies;
    double* magnitudes;
    int32_t bin_count;
} FftResult;

FFI_PLUGIN_EXPORT void generate_hanning_window(double* window, int32_t size);

FFI_PLUGIN_EXPORT FftResult* compute_real_fft(
    const double* input,
    int32_t input_length,
    double sampling_frequency
);

FFI_PLUGIN_EXPORT FftResult* compute_magnitude_spectrum(
    const double* pre_processed_input,
    int32_t input_length,
    double sampling_frequency
);

FFI_PLUGIN_EXPORT void free_fft_result(FftResult* result);

/* ── Welch's Method (welch_method.c) ── */

typedef struct {
    double* frequencies;
    double* magnitudes;
    int32_t bin_count;
} WelchResult;

FFI_PLUGIN_EXPORT WelchResult* compute_welch_spectrum(
    const double* input,
    int32_t input_length,
    double sampling_frequency,
    double frequency_resolution,
    double overlap_percent,
    int32_t weighting_type,
    double rms_range_multiplier,
    int32_t is_acc,
    double db_ref
);

FFI_PLUGIN_EXPORT void free_welch_result(WelchResult* result);

/* ── STFT Pipeline (stft_pipeline.c) ── */

typedef struct {
    double* spectrogram;
    int32_t spec_rows;
    int32_t spec_cols;
    uint8_t* pixels;
    int32_t pixel_width;
    int32_t pixel_height;
    double min_level;
    double max_level;
} StftResult;

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
);

typedef struct {
    double* filtered_data;
    int32_t filtered_rows;
    int32_t filtered_cols;
    uint8_t* pixels;
    int32_t pixel_width;
    int32_t pixel_height;
} FilteredPixelResult;

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
);

FFI_PLUGIN_EXPORT void free_filtered_pixel_result(FilteredPixelResult* result);

typedef struct {
    uint8_t* pixels;
    int32_t pixel_width;
    int32_t pixel_height;
} PixelBufferResult;

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
);

FFI_PLUGIN_EXPORT void free_pixel_buffer_result(PixelBufferResult* result);

FFI_PLUGIN_EXPORT void free_stft_result(StftResult* result);

#ifdef __cplusplus
}
#endif

#endif /* STFT_FFT_PLUGIN_H */
