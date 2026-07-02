# STFT / FFT amplitude · dB convention

이 문서는 `flutter_stft_fft_plugin` 가 산출하는 amplitude · dB 의 정확한 정의와
단위 가정을 정리한다. 호출자(앱·다른 라이브러리)가 결과 dB 를 ISO 1683 표준 또는
다른 측정 시스템의 출력과 비교할 수 있도록 한다.

---

## 1. 표준 dB 변환식

`compute_welch_spectrum` (`welch_method.c`) 과 `compute_stft_spectrogram`
(`stft_pipeline.c`) 은 동일한 dB 변환 식을 사용한다:

```
db_offset = is_acc ? 20·log10(9.807 / db_ref) : 20·log10(1.0 / db_ref)
mag_db    = db_offset + 20·log10(amplitude)
          = 20·log10( (is_acc ? 9.807 : 1.0) · amplitude / db_ref )
```

| 모드 | 풀린 식 |
|---|---|
| `is_acc = true` | `20·log10( 9.807 · amplitude / db_ref )` |
| `is_acc = false` | `20·log10( amplitude / db_ref )` |

### 단위 가정 — `is_acc = true`

`is_acc = true` 일 때 식 내부에서 `9.807` (1 g = 9.807 m/s²) 이 자동으로 곱해진다.
이것이 의미 있게 동작하려면 다음 규약을 따라야 한다:

| 인자 | 단위 |
|---|---|
| 입력 신호 (`input`) | **g** (g 단위, 1 g ≈ 9.807 m/s²) |
| `db_ref` | **m/s²** (예: 표준 진동 가속도 reference `1e-6 m/s²`) |

호출자가 다음과 같이 보내면 ISO 1683 표준 dB 와 일치하는 출력이 나온다:

```c
compute_welch_spectrum(
    /* input  */ samples_in_g,
    /* ... */
    /* is_acc */ 1,
    /* db_ref */ 1e-6      // m/s² 기준 reference, g 로 미리 환산하지 않음
);
```

⚠️ `db_ref` 를 g 단위로 미리 환산 (`1e-6 / 9.81`) 해서 보내면 `9.807` 곱셈이
이중으로 적용되어 약 +19.83 dB 만큼 오버 계산된다. **`db_ref` 는 항상 m/s² 기준
원본 reference 를 그대로 전달**한다.

### 단위 가정 — `is_acc = false`

음압(사운드) 등 다른 도메인. 입력 신호와 `db_ref` 의 단위만 일치하면 된다
(예: 입력=Pa, `db_ref = 2e-5 Pa`).

---

## 2. amplitude 정의 — RMS convention

dB 변환 식의 `amplitude` 는 **선형 RMS amplitude** 이다. 정현파 `A·cos(ω₀k)`
입력에 대해 main bin amplitude 가 `A/√2` (= `A_RMS`) 가 되도록 정규화된다.

### STFT (`stft_pipeline.c`)

```
amp_complex = sqrt(re² + im²)                         /* |X_w| */
norm        = (i != 0 && i != num_bins-1)
                ? amp_complex · (2/N)                  /* single-sided peak */
                : amp_complex · (1/N)
amplitude   = norm · √2 · weight_factor[i]            /* peak → RMS */
              ↓ apply_rms_smoothing(...)
```

Hanning window 의 coherent gain `sum(w)/N ≈ 0.5` 과 `√2` 곱셈이 결합되어
정확히 `A_RMS` 가 산출된다:
- `norm = (2/N) · |X_w| = (2/N) · A·N/4 = A/2`
- `amplitude = (A/2) · √2 = A/√2 = A_RMS` ✓

### Welch (`welch_method.c`)

```
power_sum[i] += |X_w(i)|²                              /* frame loop accumulate */
mean_pwr      = power_sum[i] / frame_count             /* time average */
factor        = (i != 0 && i != num_bins-1) ? 2.0 : 1.0  /* single-sided */
amplitude     = sqrt( mean_pwr · factor / sum(w)² )
                ↓ apply_rms_smoothing(...)
```

coherent gain `sum(w)` 기준의 정규화:
- `|X_w|² = A²·sum(w)²/4` (정현파 단일 톤)
- `mean_pwr · 2 / sum(w)² = A²/2`
- `amplitude = √(A²/2) = A/√2 = A_RMS` ✓

→ **STFT 와 Welch 는 동일 정현파 입력에 대해 동일한 amplitude (A/√2) 를 산출한다.**

### 의도된 효과

- **dB 단위 통일**: STFT 와 Welch 가 같은 신호에 대해 같은 dB 를 출력 (이전엔
  Welch 가 약 `√(N/3)` 만큼 부풀려져 dB 가 어긋났음 — §6 참조).
- **chunk_size 무관**: amplitude 가 `N` 에 의존하지 않는다. chunk_size 를 바꿔도
  같은 신호의 amplitude 와 dB 가 변하지 않는다.
- **ISO 1683 정합**: `is_acc=true` + `db_ref=1e-6` + g 단위 신호 → 표준 식
  `20·log10(9.81·A_g_RMS / 1e-6)` 과 일치.

---

## 3. ENBW 보정과 helper 식

`amplitude` 자체에는 **Hanning ENBW (= 1.5 bins) 보정이 적용되어 있지 않다**.
이는 의도된 설계 — 사용자가 좁은 대역 RMS 평균을 계산할 때 (helper 함수 또는
앱 측 후처리) ENBW 보정 인자 1.5 를 다음 식에서 한 번만 적용해야 한다:

```
df·N 범위 RMS dB = 10·log10( (Σ amplitude[i]²) / 1.5 / N_range / db_ref² )
                  ± is_acc 보정 (= 20·log10(9.807))
```

자세한 helper 식은 앱 측 `df_n_rms_helper.dart` 의 `computeDfNRmsDb` 함수와
선행 PR #80 의 ADR-003 을 참고.

### dbRef 의 단위 (helper 호출 시)

helper 식 자체는 `is_acc=true` 의 9.807 자동 곱셈을 적용하지 않는다. 따라서
호출자는 ref 측에서 단위 보정을 적용해야 한다 (앱 측 `accDbRefInG()` 패턴):

| 도착지 | 신호 단위 | 호출 시 `db_ref` |
|---|---|---|
| native (`compute_welch_spectrum`) | g | `1e-6` (raw m/s² ref) |
| native (`compute_welch_spectrum`) | m/s² | `1e-6` |
| helper (`computeDfNRmsDb`, 9.807 없음) | g | `1e-6 / 9.81` (g 단위 ref) |
| helper (9.807 없음) | m/s² | `1e-6` |

---

## 4. 호출자 사이드 sensitivity weighting

가속도계의 raw 측정값(예: `.h5` 의 `SM_Raw_Time_Data`) 이 V 단위로 저장된 경우,
native 호출 전 V → g 환산을 호출자가 적용해야 한다:

```dart
final weighting = 1000.0 / channel.sensitivity;     // 100 mV/g → 10
final samples_in_g = raw_volts.map((v) => v * weighting).toList();
```

native 는 입력 신호의 단위를 검증하지 않으며, V 가 들어오면 그대로 V 단위로 계산한다.
`is_acc=true` 의 9.807 곱셈은 *단위 환산이 아닌 g→m/s² 환산* 이라는 점을 다시 강조.

---

## 5. 검증

다음 단위 테스트가 RMS convention 을 회귀 가드한다 (`test/test_welch_method.c`,
`test/test_stft_pipeline.c`):

| 테스트 | 검증 내용 |
|---|---|
| `test_welch_rms_convention` | 정현파 A=1 입력 → amplitude main bin = 0.7071 (= 1/√2) ±5% |
| `test_welch_n_independence` | chunk_size 1024 / 2048 / 4096 에서 amplitude 동일 (RMS convention 의 N 무관성) |
| `test_welch_stft_unit_parity` | Welch peak ≈ STFT peak (±10%, 단위 통일) |
| `test_welch_amplitude_field` | amplitude → dB 재계산 invariant, byte-exact 회귀 가드 |
| `test_stft_amplitude_field` | STFT amplitude grid 가 `num_frames × num_bins` raw shape |

### 실측 검증 (앱 통합)

`100 mV/g` 가속도계 + 1 g cal 신호 (`SM_Raw_Time_Data`):

| 항목 | 이전 식 (sum(w²)) | 새 식 (RMS convention) |
|---|---|---|
| FFT Peak @ 159 Hz | 179.9 dB (V→g 환산 후) | **136.5 dB** |
| 정식 빌드 TDMS reference | — | 137.97 dB |
| 격차 | +41.93 dB | **-1.47 dB** (peak/RMS convention 미세 차이) |

→ 새 식이 표준 reference 와 정합. N=32768 가정 정확히 확정.

---

## 6. 변경 이력 — 이전 Welch amplitude 식

`v0.0.2` 이전의 Welch amplitude 는 `√(2·|X|²/sum(w²))` 식을 사용했다. 이는
정현파 A 입력 시 `A·√(N/3)` 을 산출하여 다음 문제가 있었다:

- STFT amplitude `A/√2` 와 약 `√(2N/3)` 배 차이 (N=32768 → 약 **+43.4 dB**)
- `chunk_size` 에 의존 — 같은 신호의 dB 가 frequency resolution 에 따라 변동
- ISO 1683 표준 식 `20·log10(A_RMS·9.81 / 1e-6)` 대비 N 의존성으로 오버 계산

`v0.0.2` 에서 `√(2·|X|²/sum(w)²)` (coherent gain 기준) 으로 정정. 이전 식과의
dB 격차는 `10·log10(sum(w)²/sum(w²)) = 10·log10(N · CG²/N · 1/(0.375N))` ≈
`10·log10(2N/3)` 이다 (Hanning).

### 호환성 (BREAKING)

`v0.0.2` 의 Welch `magnitudes` (dB) 와 `amplitude` (선형) 출력이 이전 버전 대비
약 `-10·log10(2N/3)` dB / `1/√(2N/3)` 배 변동한다. 호출자가 dB 임계값 또는 저장된
보고서를 사용하는 경우 마이그레이션 필요. STFT 출력은 변동 없음.

---

## 참고

- `src/welch_method.c:50-58, 100-113`
- `src/stft_pipeline.c:182-196`
- `src/signal_helpers.c:apply_rms_smoothing`
- `test/test_welch_method.c:test_welch_rms_convention` 외 2건
- 앱 측 helper: `lib/common/util/df_n_rms_helper.dart` (PR #80 / ADR-003)
- ISO 1683:2008 — Reference values for levels used in acoustics and vibration
