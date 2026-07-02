# Test Results — v0.0.2 (Welch RMS convention 정정)

`welch_method.c` 의 amplitude normalization 을 coherent gain 기준 RMS convention
으로 정정한 후의 회귀/신규 테스트 실행 결과 보존본.

- 일자: 2026-05-28
- 빌드: `cc -o <runner> <test>.c ../src/*.c ../src/kissfft/*.c -I../src -I../src/kissfft -lm -Dkiss_fft_scalar=double`
- 플랫폼: macOS (Darwin 25.2.0)

---

## test_welch_method (8/8 passed)

```
=== Welch's Method Tests ===

test_welch_basic...
  PASS (peak at 100.0 Hz, 91.0 dB)
test_welch_a_weighting...
  PASS (linear_diff=0.0, weighted_diff=30.3)
test_welch_rms_smoothing...
  PASS
test_welch_acc_db_ref...
  PASS
test_welch_amplitude_field...
  PASS
test_welch_rms_convention...
  PASS (peak_amp=0.7071, expected_rms=0.7071, ratio=1.000)
test_welch_n_independence...
  PASS (chunk 1024/2048/4096 amp = 0.7071 / 0.7071 / 0.7071, expected 0.7071)
test_welch_stft_unit_parity...
  PASS (welch_peak=0.7071, stft_peak=0.7064, ratio=1.001)

=== Results: 8 passed, 0 failed ===
```

### 핵심 검증 포인트

| 테스트 | 결과 | 해석 |
|---|---|---|
| `test_welch_rms_convention` | `peak_amp=0.7071`, `ratio=1.000` | 정현파 A=1 → amplitude = 1/√2 정확. RMS convention 확정. |
| `test_welch_n_independence` | chunk 1024/2048/4096 모두 0.7071 | amplitude 가 N 에 무관 (이전 식은 √(N/3) 의존이라 chunk 4096 에서 1.41 배). |
| `test_welch_stft_unit_parity` | Welch 0.7071 vs STFT 0.7064, ratio 1.001 | STFT 와 Welch 가 같은 단위 (±0.1% 일치). |

---

## test_stft_pipeline (30/30 passed)

```
=== STFT Pipeline Tests ===

[test_compute_stft_basic]
  PASS: result not null
  PASS: spec_rows == target_width
  PASS: spec_cols == target_height
  PASS: pixel_width > 0
  PASS: pixel_height > 0
  PASS: pixels not null
  PASS: spectrogram not null
  PASS: min < max (auto calc)

[test_filter_and_generate_pixels]
  PASS: result not null
  PASS: filtered_rows > 0
  PASS: filtered_cols > 0
  PASS: pixels not null
  PASS: pixel_width == filtered_rows
  PASS: pixel_height == filtered_cols

[test_generate_pixel_buffer]
  PASS: result not null
  PASS: pixel_width == rows
  PASS: pixel_height == cols
  PASS: pixels not null
  PASS: alpha channel is 0xFF

[test_stft_frequency_filtering]
  PASS: both results not null
  PASS: filtered height <= full height

[test_stft_amplitude_field]
  PASS: result not null
  PASS: amplitude field must not be NULL
  PASS: amp_rows must equal num_frames (15)
  PASS: amp_cols must equal num_bins (129)
  PASS: spectrogram is still resampled to target
  PASS: amplitude >= 0 and dB recomputation is finite
  PASS: baseline result not null
  PASS: baseline dimensions match result dimensions
  PASS: spectrogram must be byte-exact with baseline (regression guard)

=== Results: 30/30 passed ===
```

→ STFT 측 변경 없음. byte-exact 회귀 가드 통과.

---

## test_signal_helpers (22/22 passed)

```
PASS: linear weighting @ 1 kHz (value 0.000000000)
PASS: A-weighting @ 1 kHz (got 0.000142, tol 0.100000)
PASS: B-weighting @ 1 kHz (got 0.000339, tol 0.200000)
PASS: C-weighting @ 1 kHz (got -0.001904, tol 0.150000)
PASS: A-weighting @ 100 Hz (got -19.144954, tol 0.500000)
PASS: A-weighting @ 0 Hz guard (value -200.000000000)
PASS: smooth=1 no-op idx=0 (value 1.000000000)
PASS: smooth=1 no-op idx=1 (value 2.000000000)
PASS: smooth=1 no-op idx=2 (value 3.000000000)
PASS: smooth=1 no-op idx=3 (value 4.000000000)
PASS: smooth=1 no-op idx=4 (value 5.000000000)
PASS: smooth=0 no-op idx=0 (value 1.000000000)
PASS: smooth=0 no-op idx=1 (value 2.000000000)
PASS: smooth=0 no-op idx=2 (value 3.000000000)
PASS: smooth=0 no-op idx=3 (value 4.000000000)
PASS: smooth=0 no-op idx=4 (value 5.000000000)
PASS: smooth=3 idx=0 (got 1.581139, tol 0.000000)
PASS: smooth=3 idx=1 (got 2.160247, tol 0.000000)
PASS: smooth=3 idx=2 (got 3.109126, tol 0.000000)
PASS: smooth=3 idx=3 (got 4.082483, tol 0.000000)
PASS: smooth=3 idx=4 (got 4.527693, tol 0.000000)
PASS: smoothing empty input no-crash (value 0.000000000)

All signal_helpers tests passed.
```

→ A/B/C weighting 과 sliding RMS smoothing 동작 변화 없음.

---

## 종합

- **신규 테스트 3개** (`test_welch_rms_convention`, `test_welch_n_independence`,
  `test_welch_stft_unit_parity`) 모두 통과 — RMS convention 확정.
- **기존 회귀 가드 모두 통과** — 기존 동작 깨지지 않음 (byte-exact 비교 포함).
- 총 **60/60** 통과.

테스트 코드 위치: `test/test_welch_method.c`, `test/test_stft_pipeline.c`,
`test/test_signal_helpers.c`.

재실행 방법은 `README.md` 의 "C 테스트 실행" 섹션 참조.
