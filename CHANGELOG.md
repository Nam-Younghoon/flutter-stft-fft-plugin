## 0.0.2

### Changed (BREAKING — semantic)

* **Welch amplitude normalization 을 coherent gain (Σw) 기준 RMS convention 으로 정정.**
  이전: `amplitude = √(2·|X|²/Σw²)` (정현파 A 입력 → `A·√(N/3)`).
  이후: `amplitude = √(2·|X|²/Σw²)` (coherent gain 기준, 정현파 A 입력 → `A/√2 = A_RMS`).
  STFT (`stft_pipeline.c`) 와 동일 단위로 통일.
* 결과 `magnitudes` (dB) 와 `amplitude` (선형) 출력이 약 `-10·log10(2N/3)` dB / `1/√(2N/3)` 배
  변동. 예: `N=32768` (chunk_size) 에서 약 **-43.4 dB**. ISO 1683 표준 식과 정합.
* dB offset 식 (`is_acc` 분기, 9.807 곱셈) 은 변경 없음. STFT 측 변경 없음.

### Added

* `test/test_welch_method.c` 에 RMS convention 회귀 가드 테스트 3개 추가:
  `test_welch_rms_convention`, `test_welch_n_independence`, `test_welch_stft_unit_parity`.
* `docs/STFT_FFT_AMPLITUDE_CONVENTION.md` — amplitude / dB 정의와 단위 가정 문서화.
* `docs/test_results_v0.0.2.md` — 회귀/신규 테스트 60/60 통과 로그.

### Migration

호출자가 dB 임계값 또는 저장된 보고서를 사용하는 경우 약 `10·log10(2N/3)` dB 만큼
값이 작아진다. 새 출력은 ISO 1683 표준과 정합하므로 임계값을 재조정하거나 절대값
대신 trend 기반 비교로 전환을 권장.

## 0.0.1

* Initial scaffolding with KissFFT (double precision).
