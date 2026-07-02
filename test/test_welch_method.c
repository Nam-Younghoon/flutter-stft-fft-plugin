#define kiss_fft_scalar double
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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

/* TDD: amplitude 필드 검증 — WelchResult 에 amplitude/amp_count 추가 후 통과 예정 */
void test_welch_amplitude_field(void) {
    printf("test_welch_amplitude_field...\n");

    /* dB 변환에 사용된 상수 (welch_method.c 와 동일해야 함) */
    double ln10      = log(10.0);
    double factor    = 20.0 / ln10;
    int    is_acc    = 0;
    double db_ref    = 0.00002;
    /* is_acc=0 이므로 db_offset = factor * log(1.0 / db_ref) */
    double db_offset = factor * log(1.0 / db_ref);

    int    N  = 8192;
    double fs = 1024.0;
    double freq = 100.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = sin(2.0 * M_PI * freq * i / fs);
    }

    /* (d) 회귀 가드용 기준 결과 — weighting=0, smooth=1 (기본 경로) */
    WelchResult* baseline = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, is_acc, db_ref);

    WelchResult* result = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, is_acc, db_ref);

    /* (a) amplitude 배열 존재 */
    ASSERT_TRUE(result != NULL, "result must not be NULL");
    ASSERT_TRUE(result->amplitude != NULL,
        "amplitude field must not be NULL");

    /* (b) 크기 일치 */
    ASSERT_TRUE(result->amp_count == result->bin_count,
        "amp_count must equal bin_count");

    /* (c) amplitude → dB 재계산 값이 magnitudes 와 1e-6 이내 일치 */
    for (int i = 0; i < result->bin_count; i++) {
        double amp = result->amplitude[i];
        double expected_db = (amp > 0.0)
            ? db_offset + factor * log(amp)
            : -200.0;
        double actual_db = result->magnitudes[i];
        if (fabs(expected_db - actual_db) > 1e-6) {
            printf("  FAIL: amplitude->dB mismatch at bin=%d: expected=%.9f got=%.9f\n",
                   i, expected_db, actual_db);
            tests_failed++;
            free_welch_result(baseline);
            free_welch_result(result);
            free(input);
            return;
        }
    }

    /* (d) 기본 경로 magnitudes 가 기존 test_welch_basic 결과와 byte-exact 동일 (회귀 가드) */
    ASSERT_TRUE(baseline != NULL, "baseline result must not be NULL");
    ASSERT_TRUE(baseline->bin_count == result->bin_count,
        "baseline bin_count must match result bin_count");
    int byte_match = (memcmp(baseline->magnitudes, result->magnitudes,
        (size_t)result->bin_count * sizeof(double)) == 0);
    ASSERT_TRUE(byte_match,
        "magnitudes must be byte-exact with baseline (regression guard)");

    free_welch_result(baseline);
    free_welch_result(result);
    free(input);
    printf("  PASS\n");
    tests_passed++;
}

/* [이슈 #3] amplitude 가 RMS convention 인지 검증.
 * 정현파 A·sin(ω₀k) 입력 → main bin amplitude ≈ A/√2 (RMS).
 * 이전 식 (sum(w²) 기준) 은 A·√(N/3) 을 산출하여 STFT 와 단위 불일치였음.
 * 새 식 (coherent gain 기준) 은 STFT 와 동일 단위. */
void test_welch_rms_convention(void) {
    printf("test_welch_rms_convention...\n");

    int    N    = 8192;
    double fs   = 1024.0;
    double freq = 100.0;
    double A    = 1.0;                       /* peak amplitude */
    double expected_rms = A / sqrt(2.0);     /* ≈ 0.7071 */

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = A * sin(2.0 * M_PI * freq * i / fs);
    }

    WelchResult* result = compute_welch_spectrum(
        input, N, fs, 1.0, 75.0, 0, 1, 0, 2e-5);

    ASSERT_TRUE(result != NULL, "result must not be NULL");
    ASSERT_TRUE(result->amplitude != NULL, "amplitude must not be NULL");

    /* main bin (100 Hz 부근) peak amplitude */
    int    peak_bin = -1;
    double peak_amp = 0.0;
    for (int i = 1; i < result->bin_count - 1; i++) {
        if (result->amplitude[i] > peak_amp) {
            peak_amp = result->amplitude[i];
            peak_bin = i;
        }
    }

    /* main bin amplitude 가 A/√2 (RMS) 와 ±5% 이내.
     * Hanning leakage 로 정확히 A/√2 는 아니고 약간 작을 수 있음 (보통 96~99%). */
    double ratio = peak_amp / expected_rms;
    ASSERT_TRUE(ratio > 0.93 && ratio < 1.05,
        "main bin amplitude must equal A/sqrt(2) within +-5% (RMS convention)");

    free_welch_result(result);
    free(input);
    printf("  PASS (peak_amp=%.4f, expected_rms=%.4f, ratio=%.3f)\n",
           peak_amp, expected_rms, ratio);
    tests_passed++;
}

/* [이슈 #3] amplitude 가 chunk_size 에 무관해야 함.
 * 이전 식은 sqrt(N/3) 의존성이 있어 N 이 달라지면 dB 가 N 만큼 shift 됨.
 * 새 식은 정현파 amplitude 가 N 과 무관하게 A/√2 (RMS). */
void test_welch_n_independence(void) {
    printf("test_welch_n_independence...\n");

    double fs   = 1024.0;
    double freq = 100.0;
    double A    = 1.0;
    double expected_rms = A / sqrt(2.0);

    int N = 16384;
    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = A * sin(2.0 * M_PI * freq * i / fs);
    }

    /* fr=1 → chunk_size=1024, fr=0.5 → chunk_size=2048, fr=0.25 → chunk_size=4096 */
    WelchResult* r1 = compute_welch_spectrum(input, N, fs, 1.0,  75.0, 0, 1, 0, 2e-5);
    WelchResult* r2 = compute_welch_spectrum(input, N, fs, 0.5,  75.0, 0, 1, 0, 2e-5);
    WelchResult* r3 = compute_welch_spectrum(input, N, fs, 0.25, 75.0, 0, 1, 0, 2e-5);

    ASSERT_TRUE(r1 != NULL && r2 != NULL && r3 != NULL, "results must not be NULL");

    double a1 = 0, a2 = 0, a3 = 0;
    for (int i = 1; i < r1->bin_count - 1; i++) if (r1->amplitude[i] > a1) a1 = r1->amplitude[i];
    for (int i = 1; i < r2->bin_count - 1; i++) if (r2->amplitude[i] > a2) a2 = r2->amplitude[i];
    for (int i = 1; i < r3->bin_count - 1; i++) if (r3->amplitude[i] > a3) a3 = r3->amplitude[i];

    /* 세 chunk_size 의 peak amplitude 가 모두 A/√2 와 ±7% 이내.
     * 이전 식이라면 a3/a1 ≈ sqrt(4) = 2 배 (= +6dB) 였음. */
    ASSERT_TRUE(a1 / expected_rms > 0.93 && a1 / expected_rms < 1.07,
        "amplitude(chunk=1024) must equal A/sqrt(2)");
    ASSERT_TRUE(a2 / expected_rms > 0.93 && a2 / expected_rms < 1.07,
        "amplitude(chunk=2048) must equal A/sqrt(2)");
    ASSERT_TRUE(a3 / expected_rms > 0.93 && a3 / expected_rms < 1.07,
        "amplitude(chunk=4096) must equal A/sqrt(2)");

    free_welch_result(r1);
    free_welch_result(r2);
    free_welch_result(r3);
    free(input);
    printf("  PASS (chunk 1024/2048/4096 amp = %.4f / %.4f / %.4f, expected %.4f)\n",
           a1, a2, a3, expected_rms);
    tests_passed++;
}

/* [이슈 #3] STFT 와 Welch 가 같은 정현파 입력에 대해 같은 amplitude. */
void test_welch_stft_unit_parity(void) {
    printf("test_welch_stft_unit_parity...\n");

    int    N    = 8192;
    double fs   = 1024.0;
    double freq = 100.0;
    double A    = 1.0;

    double* input = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        input[i] = A * sin(2.0 * M_PI * freq * i / fs);
    }

    /* Welch: fr=1 → chunk_size=1024 */
    WelchResult* welch = compute_welch_spectrum(input, N, fs, 1.0, 50.0, 0, 1, 0, 2e-5);

    /* STFT: 동일 chunk_size=1024 */
    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) color_lut[i] = 0xFFFFFFFF;
    StftResult* stft = compute_stft_spectrogram(
        input, N, 1024, fs, 0, 2e-5,
        0, 50.0, 1.0,
        100, 64, 0.0, 500.0, 0.0, (double)N/fs, (double)N/fs,
        color_lut, 256, 1, 0.0, 0.0);

    ASSERT_TRUE(welch != NULL && stft != NULL, "results must not be NULL");
    ASSERT_TRUE(welch->amplitude != NULL && stft->amplitude != NULL,
        "amplitude arrays must not be NULL");

    /* Welch peak amplitude */
    double welch_peak = 0.0;
    for (int i = 1; i < welch->bin_count - 1; i++) {
        if (welch->amplitude[i] > welch_peak) welch_peak = welch->amplitude[i];
    }

    /* STFT peak amplitude (raw grid: num_frames × num_bins) */
    double stft_peak = 0.0;
    int total_cells = stft->amp_rows * stft->amp_cols;
    for (int i = 0; i < total_cells; i++) {
        if (stft->amplitude[i] > stft_peak) stft_peak = stft->amplitude[i];
    }

    /* 두 amplitude 가 ±10% 이내 (overlap/window leakage 미세 차이 허용) */
    double ratio = welch_peak / stft_peak;
    ASSERT_TRUE(ratio > 0.9 && ratio < 1.1,
        "Welch and STFT amplitude must match within +-10% (unit parity)");

    free_welch_result(welch);
    free_stft_result(stft);
    free(input);
    printf("  PASS (welch_peak=%.4f, stft_peak=%.4f, ratio=%.3f)\n",
           welch_peak, stft_peak, ratio);
    tests_passed++;
}

int main(void) {
    printf("=== Welch's Method Tests ===\n\n");

    test_welch_basic();
    test_welch_a_weighting();
    test_welch_rms_smoothing();
    test_welch_acc_db_ref();
    test_welch_amplitude_field();
    test_welch_rms_convention();
    test_welch_n_independence();
    test_welch_stft_unit_parity();

    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
