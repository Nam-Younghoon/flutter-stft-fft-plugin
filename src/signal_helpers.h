/*
 * signal_helpers.h — Welch's method 와 STFT pipeline 이 공유하는 신호 헬퍼.
 *
 * ADR-002 (sw2-noisescope-app feat/stft-fft-analysis-sync) 채택.
 * 두 파이프라인이 동일한 A/B/C weighting 공식과 RMS smoothing 로직을
 * 사용하도록 단일 소스로 분리한다. welch_method.c 에 static 으로 있던
 * 함수를 추출해 외부 visible 로 변경.
 */
#ifndef STFT_FFT_PLUGIN_SIGNAL_HELPERS_H
#define STFT_FFT_PLUGIN_SIGNAL_HELPERS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * IEC 61672 기반 주파수 가중치 계수 (dB).
 *
 * weighting_type:
 *   0 = linear (no weighting; 항상 0 dB 반환)
 *   1 = A-weighting
 *   2 = B-weighting
 *   3 = C-weighting
 *
 * f <= 0 인 경우 -200 dB guard 반환.
 */
double get_weighting_db(double frequency_hz, int32_t weighting_type);

/*
 * In-place sliding RMS smoothing.
 *
 * smooth_bins <= 1 일 때 no-op. 그 이상일 때:
 *   out[i] = sqrt( mean(amplitude[i-half .. i+half]^2) )
 *   half = smooth_bins / 2; 경계는 잘려나간 윈도우로 평균.
 *
 * num_bins <= 0 또는 amplitude == NULL 인 경우 no-op.
 */
void apply_rms_smoothing(double* amplitude, int32_t num_bins, int32_t smooth_bins);

#ifdef __cplusplus
}
#endif

#endif /* STFT_FFT_PLUGIN_SIGNAL_HELPERS_H */
