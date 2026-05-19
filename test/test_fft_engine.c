#define kiss_fft_scalar double
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../src/stft_fft_plugin.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { \
        printf("  FAIL: %s (line %d)\n", msg, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
    double _a = (a), _b = (b), _t = (tol); \
    if (fabs(_a - _b) > _t) { \
        printf("  FAIL: %s — expected %.6f, got %.6f (tol=%.6f, line %d)\n", msg, _b, _a, _t, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

void test_hanning_window(void) {
    printf("test_hanning_window...\n");
    int size = 1024;
    double* window = (double*)malloc(size * sizeof(double));
    generate_hanning_window(window, size);

    ASSERT_NEAR(window[0], 0.0, 1e-10, "window[0] should be ~0");
    ASSERT_NEAR(window[size-1], 0.0, 1e-10, "window[N-1] should be ~0");

    int mid = size / 2;
    ASSERT_NEAR(window[mid], 1.0, 0.001, "window[N/2] should be ~1.0");

    double sum = 0;
    for (int i = 0; i < size; i++) {
        ASSERT_TRUE(window[i] >= 0.0 && window[i] <= 1.0, "window values in [0,1]");
        sum += window[i];
    }
    ASSERT_TRUE(sum > 0, "window sum should be positive");

    free(window);
    printf("  PASS\n");
    tests_passed++;
}

void test_real_fft_sine_peak(void) {
    printf("test_real_fft_sine_peak...\n");

    int N = 1024;
    double fs = 1024.0;
    double freq = 100.0;
    double amplitude = 1.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = amplitude * sin(2.0 * M_PI * freq * i / fs);
    }

    FftResult* result = compute_real_fft(input, N, fs);
    ASSERT_TRUE(result != NULL, "result should not be NULL");
    ASSERT_TRUE(result->bin_count == N / 2 + 1, "bin_count should be N/2+1");

    int peak_bin = -1;
    double peak_mag = -1.0;
    for (int i = 1; i < result->bin_count - 1; i++) {
        if (result->magnitudes[i] > peak_mag) {
            peak_mag = result->magnitudes[i];
            peak_bin = i;
        }
    }

    double peak_freq = result->frequencies[peak_bin];
    ASSERT_NEAR(peak_freq, freq, 2.0, "peak frequency should be near 100 Hz");
    ASSERT_TRUE(peak_mag > 0.3, "peak magnitude should be significant");

    free_fft_result(result);
    free(input);
    printf("  PASS (peak at %.1f Hz, mag=%.4f)\n", peak_freq, peak_mag);
    tests_passed++;
}

void test_real_fft_dc(void) {
    printf("test_real_fft_dc...\n");

    int N = 256;
    double fs = 256.0;
    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = 1.0;
    }

    FftResult* result = compute_real_fft(input, N, fs);
    ASSERT_TRUE(result != NULL, "result should not be NULL");

    ASSERT_NEAR(result->frequencies[0], 0.0, 1e-10, "DC bin frequency should be 0");

    free_fft_result(result);
    free(input);
    printf("  PASS\n");
    tests_passed++;
}

void test_real_fft_two_tones(void) {
    printf("test_real_fft_two_tones...\n");

    int N = 2048;
    double fs = 2048.0;
    double f1 = 200.0, f2 = 500.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * f1 * i / fs) + 0.5 * sin(2.0 * M_PI * f2 * i / fs);
    }

    FftResult* result = compute_real_fft(input, N, fs);
    ASSERT_TRUE(result != NULL, "result should not be NULL");

    int bin1 = (int)(f1 / (fs / N));
    int bin2 = (int)(f2 / (fs / N));

    ASSERT_TRUE(result->magnitudes[bin1] > result->magnitudes[bin1 - 2], "f1 peak");
    ASSERT_TRUE(result->magnitudes[bin2] > result->magnitudes[bin2 - 2], "f2 peak");
    ASSERT_TRUE(result->magnitudes[bin1] > result->magnitudes[bin2], "f1 > f2 (amplitude 1.0 vs 0.5)");

    free_fft_result(result);
    free(input);
    printf("  PASS\n");
    tests_passed++;
}

int main(void) {
    printf("=== FFT Engine Tests ===\n\n");

    test_hanning_window();
    test_real_fft_sine_peak();
    test_real_fft_dc();
    test_real_fft_two_tones();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
