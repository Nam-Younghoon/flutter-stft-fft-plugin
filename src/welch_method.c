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

static int next_power_of_two(int n) {
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

/* get_weighting_db, apply_rms_smoothing 는 signal_helpers.c 로 추출됨 (ADR-002).
 * STFT 와 Welch 가 동일 정의를 공유하기 위한 단일 소스화. */

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
) {
    if (!input || input_length < 2 || sampling_frequency <= 0.0 ||
        frequency_resolution <= 0.0) return NULL;

    /* 1. chunk/hop size */
    int chunk_size = (int)(sampling_frequency / frequency_resolution);
    chunk_size = next_power_of_two(chunk_size);
    if (chunk_size > input_length) chunk_size = next_power_of_two(input_length / 2);
    if (chunk_size < 2) return NULL;

    double hop_ratio = 1.0 - overlap_percent / 100.0;
    int hop_size = (int)(chunk_size * hop_ratio + 0.5);
    if (hop_size < 1) hop_size = 1;
    if (hop_size > chunk_size) hop_size = chunk_size;

    int num_bins = chunk_size / 2 + 1;

    /* 2. Hanning window */
    double* window = (double*)malloc(chunk_size * sizeof(double));
    if (!window) return NULL;
    generate_hanning_window(window, chunk_size);

    /* coherent gain (CG = Σw[i]). 단일 정현파 amplitude 복원 normalization 의 기준.
     * 이전엔 window energy (Σw²) 를 썼는데, 그 식은 정현파 A 입력 시 amplitude 가
     * A·√(N/3) 만큼 부풀려져 STFT (RMS convention) 와 약 +43 dB (N=32768) 격차 발생.
     * 현 식은 정현파 A → A/√2 (RMS) 를 산출하여 STFT 와 단위 통일.
     * 자세한 유도는 docs/STFT_FFT_AMPLITUDE_CONVENTION.md 참고. */
    double coh_gain = 0.0;
    for (int i = 0; i < chunk_size; i++) {
        coh_gain += window[i];
    }

    /* 3. FFT setup */
    kiss_fftr_cfg cfg = kiss_fftr_alloc(chunk_size, 0, NULL, NULL);
    if (!cfg) { free(window); return NULL; }

    double* windowed = (double*)malloc(chunk_size * sizeof(double));
    kiss_fft_cpx* fft_out = (kiss_fft_cpx*)malloc(num_bins * sizeof(kiss_fft_cpx));
    double* power_sum = (double*)calloc(num_bins, sizeof(double));

    if (!windowed || !fft_out || !power_sum) {
        free(window); free(windowed); free(fft_out); free(power_sum);
        kiss_fft_free(cfg);
        return NULL;
    }

    /* 4. Frame loop: window → FFT → power accumulation */
    int frame_count = 0;
    for (int start = 0; start + chunk_size <= input_length; start += hop_size) {
        for (int i = 0; i < chunk_size; i++) {
            windowed[i] = input[start + i] * window[i];
        }
        kiss_fftr(cfg, windowed, fft_out);
        for (int i = 0; i < num_bins; i++) {
            power_sum[i] += fft_out[i].r * fft_out[i].r + fft_out[i].i * fft_out[i].i;
        }
        frame_count++;
    }

    free(window);
    free(windowed);
    free(fft_out);
    kiss_fft_free(cfg);

    if (frame_count == 0) {
        free(power_sum);
        return NULL;
    }

    /* 5. Average → amplitude spectrum */
    double* amplitude = (double*)malloc(num_bins * sizeof(double));
    if (!amplitude) { free(power_sum); return NULL; }

    double inv_frame_count = 1.0 / (double)frame_count;
    double inv_cg2 = (coh_gain > 0.0) ? (1.0 / (coh_gain * coh_gain)) : 1.0;

    for (int i = 0; i < num_bins; i++) {
        double mean_pwr = power_sum[i] * inv_frame_count;
        /* single-sided doubling: non-edge bin (DC/Nyquist 제외) 에만 ×2 */
        double factor = (i != 0 && i != num_bins - 1) ? 2.0 : 1.0;
        amplitude[i] = sqrt(mean_pwr * factor * inv_cg2);
        /* 정현파 A 입력 → main bin amplitude = A/√2 (RMS).
         * STFT (stft_pipeline.c) 와 동일 convention. */
    }
    free(power_sum);

    /* 6. A/B/C frequency weighting */
    if (weighting_type > 0) {
        double df = sampling_frequency / (double)chunk_size;
        for (int i = 0; i < num_bins; i++) {
            double freq = (double)i * df;
            double weight_db = get_weighting_db(freq, weighting_type);
            amplitude[i] *= pow(10.0, weight_db / 20.0);
        }
    }

    /* 7. RMS smoothing (signal_helpers 로 추출됨, ADR-002).
     * 동작은 추출 전과 동일: smooth_bins<=1 이면 no-op. */
    int smooth_bins = (int)(rms_range_multiplier);
    apply_rms_smoothing(amplitude, num_bins, smooth_bins);
    double* smoothed = amplitude;

    /* 8. dB conversion */
    WelchResult* result = (WelchResult*)malloc(sizeof(WelchResult));
    if (!result) { free(smoothed); return NULL; }

    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitudes = (double*)malloc(num_bins * sizeof(double));
    result->bin_count = num_bins;
    /* [이슈 #2 / ADR-003] free_welch_result 에서 안전하게 free 하기 위해 amplitude 를 NULL 로 초기화. */
    result->amplitude = NULL;
    result->amp_count = 0;

    if (!result->frequencies || !result->magnitudes) {
        free(smoothed);
        free_welch_result(result);
        return NULL;
    }

    double ln10 = log(10.0);
    double factor = 20.0 / ln10;
    double db_offset = is_acc
        ? factor * log(9.807 / db_ref)
        : factor * log(1.0 / db_ref);

    double df = sampling_frequency / (double)chunk_size;
    for (int i = 0; i < num_bins; i++) {
        result->frequencies[i] = (double)i * df;
        double x = smoothed[i];
        result->magnitudes[i] = (x > 0.0) ? db_offset + factor * log(x) : -200.0;
    }

    /* [이슈 #2] smoothed (= amplitude 의 별칭) 의 소유권을 result 에 양도. free 하지 않음. */
    result->amplitude = smoothed;
    result->amp_count = num_bins;
    return result;
}

FFI_PLUGIN_EXPORT void free_welch_result(WelchResult* result) {
    if (!result) return;
    free(result->frequencies);
    free(result->magnitudes);
    free(result->amplitude);
    free(result);
}
