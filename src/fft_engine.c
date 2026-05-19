#define kiss_fft_scalar double
#include "stft_fft_plugin.h"
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

FFI_PLUGIN_EXPORT void generate_hanning_window(double* window, int32_t size) {
    if (size <= 0) return;
    if (size == 1) {
        window[0] = 1.0;
        return;
    }
    double denom = (double)(size - 1);
    for (int32_t i = 0; i < size; i++) {
        window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / denom));
    }
}

FFI_PLUGIN_EXPORT FftResult* compute_real_fft(
    const double* input,
    int32_t input_length,
    double sampling_frequency
) {
    if (!input || input_length < 2 || sampling_frequency <= 0.0) return NULL;

    int nfft = input_length;
    int num_bins = nfft / 2 + 1;

    double* windowed = (double*)malloc(nfft * sizeof(double));
    if (!windowed) return NULL;

    double* window = (double*)malloc(nfft * sizeof(double));
    if (!window) { free(windowed); return NULL; }

    generate_hanning_window(window, nfft);

    for (int i = 0; i < nfft; i++) {
        windowed[i] = input[i] * window[i];
    }

    double window_energy_sum = 0.0;
    for (int i = 0; i < nfft; i++) {
        window_energy_sum += window[i] * window[i];
    }
    free(window);

    kiss_fftr_cfg cfg = kiss_fftr_alloc(nfft, 0, NULL, NULL);
    if (!cfg) { free(windowed); return NULL; }

    kiss_fft_cpx* fft_out = (kiss_fft_cpx*)malloc(num_bins * sizeof(kiss_fft_cpx));
    if (!fft_out) { kiss_fft_free(cfg); free(windowed); return NULL; }

    kiss_fftr(cfg, windowed, fft_out);
    free(windowed);
    kiss_fft_free(cfg);

    FftResult* result = (FftResult*)malloc(sizeof(FftResult));
    if (!result) { free(fft_out); return NULL; }

    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitudes = (double*)malloc(num_bins * sizeof(double));
    result->bin_count = num_bins;

    if (!result->frequencies || !result->magnitudes) {
        free_fft_result(result);
        free(fft_out);
        return NULL;
    }

    double df = sampling_frequency / (double)nfft;

    for (int i = 0; i < num_bins; i++) {
        result->frequencies[i] = (double)i * df;

        double power = (fft_out[i].r * fft_out[i].r + fft_out[i].i * fft_out[i].i);
        if (window_energy_sum > 0.0) {
            power /= window_energy_sum;
        }
        if (i != 0 && i != num_bins - 1) {
            power *= 2.0;
        }
        result->magnitudes[i] = sqrt(power);
    }

    free(fft_out);
    return result;
}

FFI_PLUGIN_EXPORT FftResult* compute_magnitude_spectrum(
    const double* pre_processed_input,
    int32_t input_length,
    double sampling_frequency
) {
    if (!pre_processed_input || input_length < 2 || sampling_frequency <= 0.0) return NULL;

    int nfft = input_length;
    int num_bins = nfft / 2 + 1;

    kiss_fftr_cfg cfg = kiss_fftr_alloc(nfft, 0, NULL, NULL);
    if (!cfg) return NULL;

    kiss_fft_cpx* fft_out = (kiss_fft_cpx*)malloc(num_bins * sizeof(kiss_fft_cpx));
    if (!fft_out) { kiss_fft_free(cfg); return NULL; }

    kiss_fftr(cfg, pre_processed_input, fft_out);
    kiss_fft_free(cfg);

    FftResult* result = (FftResult*)malloc(sizeof(FftResult));
    if (!result) { free(fft_out); return NULL; }

    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitudes = (double*)malloc(num_bins * sizeof(double));
    result->bin_count = num_bins;

    if (!result->frequencies || !result->magnitudes) {
        free_fft_result(result);
        free(fft_out);
        return NULL;
    }

    double df = sampling_frequency / (double)nfft;
    for (int i = 0; i < num_bins; i++) {
        result->frequencies[i] = (double)i * df;
        result->magnitudes[i] = sqrt(fft_out[i].r * fft_out[i].r + fft_out[i].i * fft_out[i].i);
    }

    free(fft_out);
    return result;
}

FFI_PLUGIN_EXPORT void free_fft_result(FftResult* result) {
    if (!result) return;
    free(result->frequencies);
    free(result->magnitudes);
    free(result);
}
