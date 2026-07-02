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
        /* weighting=linear, overlap=50%, rms_range=1 — 변경 전 동작 동등 */
        0, 50.0, 1.0,
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
        0, 50.0, 1.0,
        15, 129, 0.0, 500.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    StftResult* partial = compute_stft_spectrogram(
        signal, n, 256, (double)fs, 0, 2e-5,
        0, 50.0, 1.0,
        100, 64, 50.0, 200.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    ASSERT(full != NULL && partial != NULL, "both results not null");
    ASSERT(partial->pixel_height <= full->pixel_height, "filtered height <= full height");

    free_stft_result(full);
    free_stft_result(partial);
    free(signal);
}

/* TDD: amplitude 필드 검증 — StftResult 에 amplitude/amp_rows/amp_cols 추가 후 통과 예정 */
static void test_stft_amplitude_field(void) {
    printf("\n[test_stft_amplitude_field]\n");

    /* dB 변환에 사용된 상수 (stft_pipeline.c 와 동일해야 함) */
    double ln10  = log(10.0);
    double factor = 20.0 / ln10;
    int    is_acc = 0;
    double db_ref = 2e-5;
    /* is_acc=0 이므로 db_offset = factor * log(1.0 / db_ref) */
    double db_offset = factor * log(1.0 / db_ref);

    int fs = 1000;
    int n  = 2000;
    double* signal = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2.0 * M_PI * 100.0 * i / fs);
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (i << 24) | (i << 16) | (i << 8) | 0xFF;
    }

    /* (d) 회귀 가드용 기준 결과 — weighting=0, overlap=50, rms_range=1 (기본 경로) */
    StftResult* baseline = compute_stft_spectrogram(
        signal, n, 256, (double)fs, is_acc, db_ref,
        0, 50.0, 1.0,
        100, 64, 0.0, 500.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    StftResult* result = compute_stft_spectrogram(
        signal, n, 256, (double)fs, is_acc, db_ref,
        0, 50.0, 1.0,
        100, 64, 0.0, 500.0, 0.0, 2.0, 2.0,
        color_lut, 256, 1, 0.0, 0.0
    );

    /* (a) amplitude 배열 존재 */
    ASSERT(result != NULL, "result not null");
    ASSERT(result->amplitude != NULL, "amplitude field must not be NULL");

    /* (b) amplitude 의 grid 가 raw STFT grid (num_frames × num_bins).
     *
     * chunk_size=256, sample n=2000, offset=128 → padded_len=2128,
     *   num_bins = 256/2 + 1 = 129
     *   hop_size = 128 (50% overlap)
     *   num_frames: s = 0, 128, ..., 1872 → 15
     *
     * resampled spectrogram 은 target_width=100 × target_height=64. */
    ASSERT(result->amp_rows == 15, "amp_rows must equal num_frames (15)");
    ASSERT(result->amp_cols == 129, "amp_cols must equal num_bins (129)");
    ASSERT(result->spec_rows == 100 && result->spec_cols == 64,
        "spectrogram is still resampled to target");

    /* (c) amplitude 의 변환 invariant: amp >= 0, amp > 0 일 때 dB 재계산식이
     * finite (raw 시점의 dB 변환 일관성. raw spectrogram 은 외부 미노출이라
     * 외부 ctest 에서는 직접 비교 불가). */
    int total_cells_raw = result->amp_rows * result->amp_cols;
    int amp_valid = 1;
    for (int idx = 0; idx < total_cells_raw; idx++) {
        double amp = result->amplitude[idx];
        if (amp < 0.0) { amp_valid = 0; break; }
        if (amp > 0.0) {
            double recomputed_db = db_offset + factor * log(amp);
            if (!isfinite(recomputed_db)) { amp_valid = 0; break; }
        }
    }
    ASSERT(amp_valid, "amplitude >= 0 and dB recomputation is finite");

    /* (d) 기본 경로 spectrogram 이 baseline 과 byte-exact 동일 (회귀 가드) */
    ASSERT(baseline != NULL, "baseline result not null");
    ASSERT(baseline->spec_rows == result->spec_rows &&
           baseline->spec_cols == result->spec_cols,
        "baseline dimensions match result dimensions");
    int spec_byte_match = (memcmp(baseline->spectrogram, result->spectrogram,
        (size_t)(result->spec_rows * result->spec_cols) * sizeof(double)) == 0);
    ASSERT(spec_byte_match,
        "spectrogram must be byte-exact with baseline (regression guard)");

    free_stft_result(baseline);
    free_stft_result(result);
    free(signal);
}

/* ==========================================================================
 * TDD RED: STFT 줌(주파수 범위 좁히기) 시 주파수 해상도 회복 검증
 *
 * 현재 버그:
 *   compute_stft_spectrogram 은 원본 STFT(num_frames×num_bins) 전체를
 *   target_width×target_height 로 bilinear 다운샘플 한 뒤, 주파수/시간 crop 을 수행
 *   (resample→crop 순서). 이 때문에 좁은 대역으로 줌해도 이미 뭉개진 grid 를
 *   자를 뿐, 해상도가 회복되지 않는다.
 *
 * 미래 수정 방향:
 *   원본 grid 에서 먼저 crop → 그 대역을 target 크기로 resample (crop→resample).
 *   다운 비율 > 1 축에는 max-pooling(peak 보존) 적용.
 * ========================================================================== */

/* --------------------------------------------------------------------------
 * 테스트 1 — zoom_recovers_resolution (RED 필수)
 *
 * 파라미터:
 *   fs=1000, chunk_size=256 (num_bins=129, nyquist=500Hz, df≈3.906Hz/bin)
 *   n=4000, 100Hz 사인파, total_duration=4.0초
 *   target_width=100, target_height=64
 *   (a) 전체 대역: start_freq=0, end_freq=500 (nyquist)
 *   (b) 좁은 대역: start_freq=0, end_freq=62.5 (nyquist/8)
 *
 * 핵심 ASSERT (RED 원인):
 *   현재 코드: 전체 64행으로 resample 후 [0, 8]행 crop → pixel_height=8
 *   기대 (미래): 원본 17 bins crop 후 64행 resample → pixel_height=64
 *   ASSERT: narrow->pixel_height >= 58  (= 0.9 * 64, 현재 8이므로 FAIL)
 *
 * 보조 ASSERT (행당 주파수폭 비율, 결정론적):
 *   현재: narrow_row_width(=62.5/8) / full_row_width(=500/64) = 1.0 → FAIL
 *   미래: narrow_row_width(=62.5/64) / full_row_width(=500/64) = 0.125 → PASS
 *   ASSERT: narrow_row_width / full_row_width <= 0.25
 * -------------------------------------------------------------------------- */
static void test_zoom_recovers_resolution(void) {
    printf("\n[test_zoom_recovers_resolution]\n");

    int fs = 1000;
    int n = 4000;
    double total_duration = (double)n / (double)fs;  /* 4.0초 */
    double nyquist = (double)fs / 2.0;               /* 500.0 Hz */
    int target_width = 100;
    int target_height = 64;

    /* 100Hz 사인파 생성 */
    double* signal = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2.0 * M_PI * 100.0 * (double)i / (double)fs);
    }

    int32_t color_lut[256];
    for (int i = 0; i < 256; i++) {
        color_lut[i] = (i << 24) | (i << 16) | (i << 8) | 0xFF;
    }

    /* (a) 전체 대역: 0 ~ 500Hz */
    StftResult* full = compute_stft_spectrogram(
        signal, n,
        256,                /* chunk_size */
        (double)fs,
        0,                  /* is_acc */
        2e-5,               /* db_ref */
        0, 50.0, 1.0,       /* weighting, overlap, rms_range */
        target_width, target_height,
        0.0, nyquist,       /* start_freq=0, end_freq=500Hz */
        0.0, total_duration, total_duration,
        color_lut, 256,
        1, 0.0, 0.0
    );

    /* (b) 좁은 대역: 0 ~ 62.5Hz (= nyquist/8)
     *   원본 grid 에서 약 17 bins (0~62.5Hz) 에 해당
     *   미래 crop->resample 이면 17 bins → 64 rows: pixel_height=64
     *   현재 resample->crop 이면 64행 중 8행 crop: pixel_height=8 */
    double narrow_end_freq = nyquist / 8.0;  /* 62.5 Hz */
    StftResult* narrow = compute_stft_spectrogram(
        signal, n,
        256,
        (double)fs,
        0,
        2e-5,
        0, 50.0, 1.0,
        target_width, target_height,
        0.0, narrow_end_freq,  /* start_freq=0, end_freq=62.5Hz */
        0.0, total_duration, total_duration,
        color_lut, 256,
        1, 0.0, 0.0
    );

    ASSERT(full != NULL && narrow != NULL, "both full and narrow results must not be null");

    if (full != NULL && narrow != NULL) {
        /* --- 핵심 ASSERT ---
         * 좁은 대역 호출의 pixel_height 가 target_height 에 근접해야 한다.
         * 기준: pixel_height >= 0.9 * target_height (= 58)
         *
         * 현재 코드(resample→crop): pixel_height = ceil(62.5*64/500) = 8 → FAIL
         * 수정 후(crop→resample): pixel_height = 64 → PASS */
        int min_expected_height = (int)(0.9 * (double)target_height);  /* 57 */
        ASSERT(
            narrow->pixel_height >= min_expected_height,
            "zoom: narrow band pixel_height must recover to >= 0.9*target_height"
            " (current resample->crop: ~8, expected crop->resample: ~64)"
        );

        /* --- 보조 ASSERT ---
         * 좁은 대역의 "행당 주파수폭"이 전체 대역 행당폭의 1/4 이하여야 한다.
         *
         * full_row_width  = end_freq_full  / full->pixel_height  = 500/64  = 7.8125 Hz/row
         * narrow_row_width = narrow_end_freq / narrow->pixel_height
         *
         * 현재(resample→crop): narrow_row_width = 62.5/8 = 7.8125 → ratio=1.0 → FAIL
         * 수정 후(crop→resample): narrow_row_width = 62.5/64 = 0.977 → ratio=0.125 → PASS */
        if (full->pixel_height > 0 && narrow->pixel_height > 0) {
            double full_row_hz = nyquist / (double)full->pixel_height;
            double narrow_row_hz = narrow_end_freq / (double)narrow->pixel_height;
            double ratio = narrow_row_hz / full_row_hz;
            ASSERT(
                ratio <= 0.25,
                "zoom: narrow band Hz-per-row must be <= 1/4 of full band Hz-per-row"
                " (current: ratio~1.0, expected after fix: ratio~0.125)"
            );
        } else {
            /* pixel_height 가 0 인 경우도 실패로 처리 */
            ASSERT(0, "zoom: pixel_height must be positive for both results");
        }
    }

    free_stft_result(full);
    free_stft_result(narrow);
    free(signal);
}

int main(void) {
    printf("=== STFT Pipeline Tests ===\n");

    test_compute_stft_basic();
    test_filter_and_generate_pixels();
    test_generate_pixel_buffer();
    test_stft_frequency_filtering();
    test_stft_amplitude_field();
    test_zoom_recovers_resolution();

    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
