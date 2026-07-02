/*
 * signal_helpers.c — get_weighting_db / apply_rms_smoothing 구현.
 *
 * ADR-002: 본 함수들은 welch_method.c 에서 추출되어 STFT/Welch 양쪽에서
 * 공유된다. 동작은 추출 전과 byte-exact 동일해야 한다 (회귀 가드).
 */
#include "signal_helpers.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* get_weighting_db — IEC 61672 A/B/C weighting (dB).
 *
 * 원 코드: welch_method.c line 18~54 (static 함수). static 키워드를 제거하고
 * external linkage 로 옮김. 공식과 정규화 상수(2.0 / 0.17 / 0.06) 모두 동일.
 */
double get_weighting_db(double frequency_hz, int32_t weighting_type) {
    if (frequency_hz <= 0.0) return -200.0;

    double f = frequency_hz;
    double f2 = f * f;
    double f4 = f2 * f2;
    double c12194_sq = 12194.0 * 12194.0;

    switch (weighting_type) {
    case 1: { /* A-weighting */
        double num = c12194_sq * f4;
        double den = (f2 + 20.6 * 20.6) *
                     (f2 + c12194_sq) *
                     sqrt((f2 + 107.7 * 107.7) * (f2 + 737.9 * 737.9));
        if (den == 0.0) return -200.0;
        double ra = num / den;
        return 20.0 * log10(ra) + 2.0;
    }
    case 2: { /* B-weighting */
        double num = c12194_sq * f * f2;
        double den = (f2 + 20.6 * 20.6) *
                     (f2 + c12194_sq) *
                     sqrt(f2 + 158.5 * 158.5);
        if (den == 0.0) return -200.0;
        double rb = num / den;
        return 20.0 * log10(rb) + 0.17;
    }
    case 3: { /* C-weighting */
        double num = c12194_sq * f2;
        double den = (f2 + 20.6 * 20.6) * (f2 + c12194_sq);
        if (den == 0.0) return -200.0;
        double rc = num / den;
        return 20.0 * log10(rc) + 0.06;
    }
    default:
        return 0.0; /* linear */
    }
}

/* apply_rms_smoothing — in-place sliding RMS.
 *
 * 원 코드: welch_method.c line 157~178 (외부 malloc + 결과 swap 형태).
 * 추출하면서 호출자가 별도 buffer 를 관리하지 않도록 in-place 로 통일.
 * 동작 의미는 동일 — smooth<=1 일 때 입력 그대로, 그 외에는 sliding window
 * RMS. 경계는 잘려나간 윈도우 (count 만큼만 평균).
 */
void apply_rms_smoothing(double* amplitude, int32_t num_bins, int32_t smooth_bins) {
    if (amplitude == NULL || num_bins <= 0 || smooth_bins <= 1) return;

    double* tmp = (double*)malloc((size_t)num_bins * sizeof(double));
    if (tmp == NULL) return;

    int half_bins = smooth_bins / 2;
    for (int i = 0; i < num_bins; i++) {
        double sum_sq = 0.0;
        int count = 0;
        for (int j = i - half_bins; j <= i + half_bins; j++) {
            if (j >= 0 && j < num_bins) {
                sum_sq += amplitude[j] * amplitude[j];
                count++;
            }
        }
        tmp[i] = (count > 0) ? sqrt(sum_sq / (double)count) : 0.0;
    }
    memcpy(amplitude, tmp, (size_t)num_bins * sizeof(double));
    free(tmp);
}
