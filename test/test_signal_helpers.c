/*
 * test_signal_helpers.c — signal_helpers 의 단위 테스트
 *
 * 검증 대상:
 *   1. get_weighting_db(f, type) — IEC 61672 A/B/C-weighting reference (1 kHz @ 0 dB)
 *   2. apply_rms_smoothing(amp, n, smooth_bins) — in-place sliding RMS, smooth<=1 no-op
 *
 * ADR-002 의 헬퍼 추출 검증. welch_method.c 의 동작과 동등해야 한다.
 */

#include "signal_helpers.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define EXPECT_NEAR(actual, expected, tol, msg) \
    do { \
        double a = (actual); \
        double e = (expected); \
        if (fabs(a - e) > (tol)) { \
            fprintf(stderr, "FAIL: %s — expected %.6f, got %.6f (tol %.6f)\n", \
                    msg, e, a, (double)(tol)); \
            failures++; \
        } else { \
            fprintf(stdout, "PASS: %s (got %.6f, tol %.6f)\n", msg, a, (double)(tol)); \
        } \
    } while (0)

#define EXPECT_EQ_DOUBLE(actual, expected, msg) \
    do { \
        double a = (actual); \
        double e = (expected); \
        if (a != e) { \
            fprintf(stderr, "FAIL: %s — expected %.9f, got %.9f\n", msg, e, a); \
            failures++; \
        } else { \
            fprintf(stdout, "PASS: %s (value %.9f)\n", msg, a); \
        } \
    } while (0)

/* ── A/B/C weighting reference (IEC 61672, 1 kHz @ 0 dB) ─────────────────── */

static void test_weighting_reference_1kHz(void) {
    /* linear (type=0): 항상 0 dB */
    EXPECT_EQ_DOUBLE(get_weighting_db(1000.0, 0), 0.0,
                     "linear weighting @ 1 kHz");

    /* A-weighting (type=1): 1 kHz @ 0 dB ± 0.1 (IEC 61672) */
    EXPECT_NEAR(get_weighting_db(1000.0, 1), 0.0, 0.1,
                "A-weighting @ 1 kHz");

    /* B-weighting (type=2): 1 kHz @ 0.17 dB normalization */
    EXPECT_NEAR(get_weighting_db(1000.0, 2), 0.17, 0.2,
                "B-weighting @ 1 kHz");

    /* C-weighting (type=3): 1 kHz @ 0.06 dB normalization */
    EXPECT_NEAR(get_weighting_db(1000.0, 3), 0.06, 0.15,
                "C-weighting @ 1 kHz");
}

static void test_weighting_low_freq_attenuation(void) {
    /* A-weighting @ 100 Hz: 약 -19.1 dB (IEC 61672 표) */
    EXPECT_NEAR(get_weighting_db(100.0, 1), -19.1, 0.5,
                "A-weighting @ 100 Hz");

    /* A-weighting @ 0 Hz: -200 (guard) */
    EXPECT_EQ_DOUBLE(get_weighting_db(0.0, 1), -200.0,
                     "A-weighting @ 0 Hz guard");
}

/* ── apply_rms_smoothing ─────────────────────────────────────────────────── */

static void test_smoothing_noop_smooth_le_1(void) {
    double amp[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double original[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

    apply_rms_smoothing(amp, 5, 1);
    for (int i = 0; i < 5; i++) {
        char buf[64];
        snprintf(buf, sizeof buf, "smooth=1 no-op idx=%d", i);
        EXPECT_EQ_DOUBLE(amp[i], original[i], buf);
    }

    /* smooth=0 도 no-op */
    apply_rms_smoothing(amp, 5, 0);
    for (int i = 0; i < 5; i++) {
        char buf[64];
        snprintf(buf, sizeof buf, "smooth=0 no-op idx=%d", i);
        EXPECT_EQ_DOUBLE(amp[i], original[i], buf);
    }
}

static void test_smoothing_window_3(void) {
    /* sliding RMS, window=3 (half=1): out[i] = sqrt(mean(amp[i-1..i+1]^2)) */
    double amp[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    apply_rms_smoothing(amp, 5, 3);

    /* idx=0: [1, 2] (j=-1 skip) → sqrt((1+4)/2) = sqrt(2.5) ≈ 1.5811 */
    EXPECT_NEAR(amp[0], sqrt(2.5), 1e-9, "smooth=3 idx=0");
    /* idx=1: [1, 2, 3] → sqrt((1+4+9)/3) = sqrt(14/3) ≈ 2.1602 */
    EXPECT_NEAR(amp[1], sqrt(14.0 / 3.0), 1e-9, "smooth=3 idx=1");
    /* idx=2: [2, 3, 4] → sqrt((4+9+16)/3) = sqrt(29/3) ≈ 3.1091 */
    EXPECT_NEAR(amp[2], sqrt(29.0 / 3.0), 1e-9, "smooth=3 idx=2");
    /* idx=3: [3, 4, 5] → sqrt((9+16+25)/3) = sqrt(50/3) ≈ 4.0825 */
    EXPECT_NEAR(amp[3], sqrt(50.0 / 3.0), 1e-9, "smooth=3 idx=3");
    /* idx=4: [4, 5] → sqrt((16+25)/2) = sqrt(20.5) ≈ 4.5277 */
    EXPECT_NEAR(amp[4], sqrt(20.5), 1e-9, "smooth=3 idx=4");
}

static void test_smoothing_empty_input(void) {
    /* num_bins=0: no crash, no work */
    apply_rms_smoothing(NULL, 0, 3);
    EXPECT_EQ_DOUBLE(0.0, 0.0, "smoothing empty input no-crash");
}

int main(void) {
    test_weighting_reference_1kHz();
    test_weighting_low_freq_attenuation();
    test_smoothing_noop_smooth_le_1();
    test_smoothing_window_3();
    test_smoothing_empty_input();

    if (failures > 0) {
        fprintf(stderr, "\n%d test(s) failed.\n", failures);
        return 1;
    }
    fprintf(stdout, "\nAll signal_helpers tests passed.\n");
    return 0;
}
