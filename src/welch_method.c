#define kiss_fft_scalar double
#include "stft_fft_plugin.h"
#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fftr.h"
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

static double get_weighting_db(double f, int32_t weighting_type) {
    if (f <= 0.0) return -200.0;

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
        return 0.0;
    }
}

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

    double window_energy_sum = 0.0;
    for (int i = 0; i < chunk_size; i++) {
        window_energy_sum += window[i] * window[i];
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
    double inv_window_energy = (window_energy_sum > 0.0) ? (1.0 / window_energy_sum) : 1.0;

    for (int i = 0; i < num_bins; i++) {
        double power = power_sum[i] * inv_frame_count * inv_window_energy;
        if (i != 0 && i != num_bins - 1) {
            power *= 2.0;
        }
        amplitude[i] = sqrt(power);
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

    /* 7. RMS smoothing */
    int smooth_bins = (int)(rms_range_multiplier);
    double* smoothed;
    if (smooth_bins > 1) {
        smoothed = (double*)malloc(num_bins * sizeof(double));
        if (!smoothed) { free(amplitude); return NULL; }

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
            smoothed[i] = (count > 0) ? sqrt(sum_sq / (double)count) : 0.0;
        }
        free(amplitude);
    } else {
        smoothed = amplitude;
    }

    /* 8. dB conversion */
    WelchResult* result = (WelchResult*)malloc(sizeof(WelchResult));
    if (!result) { free(smoothed); return NULL; }

    result->frequencies = (double*)malloc(num_bins * sizeof(double));
    result->magnitudes = (double*)malloc(num_bins * sizeof(double));
    result->bin_count = num_bins;

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

    free(smoothed);
    return result;
}

FFI_PLUGIN_EXPORT void free_welch_result(WelchResult* result) {
    if (!result) return;
    free(result->frequencies);
    free(result->magnitudes);
    free(result);
}
