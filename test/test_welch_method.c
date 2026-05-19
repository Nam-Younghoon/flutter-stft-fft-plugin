#define kiss_fft_scalar double
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

void test_welch_basic(void) {
    printf("test_welch_basic...\n");

    int N = 8192;
    double fs = 1024.0;
    double freq = 100.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * freq * i / fs);
    }

    WelchResult* result = compute_welch_spectrum(
        input, N, fs,
        1.0,   /* frequency_resolution = 1 Hz */
        75.0,  /* overlap 75% */
        0,     /* weighting: linear */
        1,     /* rms_range_multiplier = 1 (no smoothing) */
        0,     /* is_acc = false */
        0.00002 /* db_ref */
    );

    ASSERT_TRUE(result != NULL, "result should not be NULL");
    ASSERT_TRUE(result->bin_count > 0, "bin_count should be positive");

    int peak_bin = -1;
    double peak_mag = -9999.0;
    for (int i = 1; i < result->bin_count - 1; i++) {
        if (result->magnitudes[i] > peak_mag) {
            peak_mag = result->magnitudes[i];
            peak_bin = i;
        }
    }

    double peak_freq = result->frequencies[peak_bin];
    ASSERT_NEAR(peak_freq, freq, 2.0, "peak frequency should be ~100 Hz");
    ASSERT_TRUE(peak_mag > 80.0, "peak dB should be > 80 (loud sine at 0dBFS)");

    free_welch_result(result);
    free(input);
    printf("  PASS (peak at %.1f Hz, %.1f dB)\n", peak_freq, peak_mag);
    tests_passed++;
}

void test_welch_a_weighting(void) {
    printf("test_welch_a_weighting...\n");

    int N = 8192;
    double fs = 8192.0;
    double freq_low = 50.0;
    double freq_high = 1000.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * freq_low * i / fs) +
                   sin(2.0 * M_PI * freq_high * i / fs);
    }

    WelchResult* linear = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, 0, 0.00002);
    WelchResult* a_weighted = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 1, 1, 0, 0.00002);

    ASSERT_TRUE(linear != NULL, "linear result should not be NULL");
    ASSERT_TRUE(a_weighted != NULL, "a_weighted result should not be NULL");

    int bin_low = (int)(freq_low);
    int bin_high = (int)(freq_high);

    double linear_diff = linear->magnitudes[bin_high] - linear->magnitudes[bin_low];
    double weighted_diff = a_weighted->magnitudes[bin_high] - a_weighted->magnitudes[bin_low];

    ASSERT_TRUE(weighted_diff > linear_diff,
        "A-weighting should boost 1kHz relative to 50Hz more than linear");

    free_welch_result(linear);
    free_welch_result(a_weighted);
    free(input);
    printf("  PASS (linear_diff=%.1f, weighted_diff=%.1f)\n", linear_diff, weighted_diff);
    tests_passed++;
}

void test_welch_rms_smoothing(void) {
    printf("test_welch_rms_smoothing...\n");

    int N = 8192;
    double fs = 1024.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * 100.0 * i / fs) +
                   0.1 * sin(2.0 * M_PI * 102.0 * i / fs);
    }

    WelchResult* no_smooth = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, 0, 0.00002);
    WelchResult* smoothed = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 5, 0, 0.00002);

    ASSERT_TRUE(no_smooth != NULL && smoothed != NULL, "results should not be NULL");

    double no_smooth_var = 0.0, smoothed_var = 0.0;
    for (int i = 10; i < no_smooth->bin_count - 10; i++) {
        double d1 = no_smooth->magnitudes[i] - no_smooth->magnitudes[i-1];
        double d2 = smoothed->magnitudes[i] - smoothed->magnitudes[i-1];
        no_smooth_var += d1 * d1;
        smoothed_var += d2 * d2;
    }

    ASSERT_TRUE(smoothed_var <= no_smooth_var,
        "smoothed spectrum should be less variable than unsmoothed");

    free_welch_result(no_smooth);
    free_welch_result(smoothed);
    free(input);
    printf("  PASS\n");
    tests_passed++;
}

void test_welch_acc_db_ref(void) {
    printf("test_welch_acc_db_ref...\n");

    int N = 4096;
    double fs = 1024.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * 100.0 * i / fs);
    }

    WelchResult* mic = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, 0, 0.00002);
    WelchResult* acc = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, 1, 0.00002);

    ASSERT_TRUE(mic != NULL && acc != NULL, "results should not be NULL");

    ASSERT_TRUE(acc->magnitudes[100] != mic->magnitudes[100],
        "acc and mic should have different dB offsets");

    free_welch_result(mic);
    free_welch_result(acc);
    free(input);
    printf("  PASS\n");
    tests_passed++;
}

int main(void) {
    printf("=== Welch's Method Tests ===\n\n");

    test_welch_basic();
    test_welch_a_weighting();
    test_welch_rms_smoothing();
    test_welch_acc_db_ref();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
