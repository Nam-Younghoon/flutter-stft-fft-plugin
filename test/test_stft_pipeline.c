#define kiss_fft_scalar double
#include "../src/stft_fft_plugin.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define ASSERT(cond, msg) do { \
    tests_run++; \
    if (!(cond)) { printf("  FAIL: %s\n", msg); } \
    else { tests_passed++; printf("  PASS: %s\n", msg); } \
} while(0)

static void test_compute_stft_basic(void) {
    printf("\n[test_compute_stft_basic]\n");

    int fs = 1000;
    int n = 2000;
    double* signal = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2.0 * M_PI * 100.0 * i / fs);
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (i << 24) | (i << 16) | (i << 8) | 0xFF;
    }

    StftResult* result = compute_stft_spectrogram(
        signal, n,
        256,
        (double)fs,
        0,
        2e-5,
        100, 64,
        0.0, 500.0,
        0.0, 2.0,
        2.0,
        color_lut, 256,
        1,
        0.0, 0.0
    );

    ASSERT(result != NULL, "result not null");
    ASSERT(result->spec_rows == 100, "spec_rows == target_width");
    ASSERT(result->spec_cols == 64, "spec_cols == target_height");
    ASSERT(result->pixel_width > 0, "pixel_width > 0");
    ASSERT(result->pixel_height > 0, "pixel_height > 0");
    ASSERT(result->pixels != NULL, "pixels not null");
    ASSERT(result->spectrogram != NULL, "spectrogram not null");
    ASSERT(result->min_level < result->max_level, "min < max (auto calc)");

    free_stft_result(result);
    free(signal);
}

static void test_filter_and_generate_pixels(void) {
    printf("\n[test_filter_and_generate_pixels]\n");

    int rows = 50;
    int cols = 32;
    double* spec = (double*)malloc(rows * cols * sizeof(double));
    for (int i = 0; i < rows * cols; i++) {
        spec[i] = -60.0 + (double)(i % 100);
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (0xFF << 24) | (i << 16) | ((255 - i) << 8) | 0xFF;
    }

    FilteredPixelResult* result = filter_and_generate_pixels(
        spec, rows, cols,
        1000.0,
        256,
        2.0,
        0.0, 500.0,
        0.0, 2.0,
        50, 32,
        color_lut, 256,
        -60.0, 40.0
    );

    ASSERT(result != NULL, "result not null");
    ASSERT(result->filtered_rows > 0, "filtered_rows > 0");
    ASSERT(result->filtered_cols > 0, "filtered_cols > 0");
    ASSERT(result->pixels != NULL, "pixels not null");
    ASSERT(result->pixel_width == result->filtered_rows, "pixel_width == filtered_rows");
    ASSERT(result->pixel_height == result->filtered_cols, "pixel_height == filtered_cols");

    free_filtered_pixel_result(result);
    free(spec);
}

static void test_generate_pixel_buffer(void) {
    printf("\n[test_generate_pixel_buffer]\n");

    int rows = 10;
    int cols = 8;
    double* data = (double*)malloc(rows * cols * sizeof(double));
    for (int i = 0; i < rows * cols; i++) {
        data[i] = -20.0 + 0.5 * i;
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (i << 24) | (i << 16) | (i << 8) | 0xFF;
    }

    PixelBufferResult* result = generate_pixel_buffer(
        data, rows, cols,
        rows, cols,
        color_lut, 256,
        -20.0, 20.0
    );

    ASSERT(result != NULL, "result not null");
    ASSERT(result->pixel_width == rows, "pixel_width == rows");
    ASSERT(result->pixel_height == cols, "pixel_height == cols");
    ASSERT(result->pixels != NULL, "pixels not null");

    int pixel_idx = (rows * (cols - 1 - 0) + 0) * 4;
    ASSERT(result->pixels[pixel_idx + 3] == 0xFF, "alpha channel is 0xFF");

    free_pixel_buffer_result(result);
    free(data);
}

static void test_stft_frequency_filtering(void) {
    printf("\n[test_stft_frequency_filtering]\n");

    int fs = 1000;
    int n = 2000;
    double* signal = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2.0 * M_PI * 100.0 * i / fs);
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (i << 24) | (i << 16) | (i << 8) | 0xFF;
    }

    StftResult* full = compute_stft_spectrogram(
        signal, n, 256, (double)fs, 0, 2e-5,
        100, 64, 0.0, 500.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    StftResult* partial = compute_stft_spectrogram(
        signal, n, 256, (double)fs, 0, 2e-5,
        100, 64, 50.0, 200.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    ASSERT(full != NULL && partial != NULL, "both results not null");
    ASSERT(partial->pixel_height <= full->pixel_height, "filtered height <= full height");

    free_stft_result(full);
    free_stft_result(partial);
    free(signal);
}

int main(void) {
    printf("=== STFT Pipeline Tests ===\n");

    test_compute_stft_basic();
    test_filter_and_generate_pixels();
    test_generate_pixel_buffer();
    test_stft_frequency_filtering();

    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
