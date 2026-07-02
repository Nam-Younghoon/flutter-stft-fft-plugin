# flutter_stft_fft_plugin

**한국어** | [English](README.en.md)

KissFFT 기반의 네이티브 FFT/STFT Flutter 플러그인입니다. Dart FFI를 통해 C로 구현된 고성능 신호처리 알고리즘을 직접 호출합니다.

## 주요 특징

- **배정밀도(Double Precision) FFT** — KissFFT 라이브러리 기반의 정밀한 주파수 분석
- **Welch's Method** — 중첩 윈도우 평균화를 통한 안정적인 전력 스펙트럼 추정
- **STFT 스펙트로그램** — 시간-주파수 2D 분석 및 실시간 픽셀 렌더링
- **A/B/C 주파수 가중치** — 음향 측정 국제 표준 지원
- **Dart FFI 직접 호출** — 플랫폼 채널 오버헤드 없는 네이티브 성능

## 지원 플랫폼

| 플랫폼 | 최소 버전 | 빌드 시스템 |
|--------|----------|------------|
| iOS | 12.0 | CocoaPods |
| Android | API 24 | CMake |

## 참고 문서

- [STFT / FFT amplitude · dB convention](docs/STFT_FFT_AMPLITUDE_CONVENTION.md) — amplitude 정의, dB 변환식, `isAcc` 단위 가정, ENBW 보정 위치 등 결과 dB 를 표준(ISO 1683)과 비교할 때 알아야 할 규약.
- [v0.0.2 테스트 결과](docs/test_results_v0.0.2.md) — Welch RMS convention 정정 후 회귀/신규 테스트 통과 로그.

## 설치

`pubspec.yaml`에 의존성을 추가합니다:

```yaml
dependencies:
  flutter_stft_fft_plugin:
    git:
      url: https://github.com/Nam-Younghoon/flutter-stft-fft-plugin.git
      ref: main
```

## API

모든 메서드는 `StftFftPlugin` 클래스의 static 메서드로 제공됩니다.

### Real FFT

시간 도메인 신호를 주파수 도메인으로 변환합니다. Hanning 윈도우가 자동 적용됩니다.

```dart
import 'package:flutter_stft_fft_plugin/stft_fft_plugin.dart';

final result = StftFftPlugin.computeRealFft(
  samples: audioSamples,        // Float64List
  samplingFrequency: 44100.0,
);

final frequencies = result['frequencies']!; // Float64List — 주파수 축 (Hz)
final magnitudes = result['magnitudes']!;   // Float64List — 진폭
```

### Welch's Method

중첩 윈도우 기반 평균화로 노이즈를 줄인 전력 스펙트럼을 계산합니다. 결과는 dB 스케일로 반환됩니다.

```dart
final result = StftFftPlugin.computeWelchSpectrum(
  samples: audioSamples,
  samplingFrequency: 44100.0,
  frequencyResolution: 1.0,    // Hz 단위 주파수 분해능
  overlapPercent: 50.0,        // 윈도우 중첩률 (%)
  weightingType: 1,            // 0=없음, 1=A, 2=B, 3=C 가중치
  rmsRangeMultiplier: 1.0,     // RMS 스무싱 범위
  isAcc: false,                // true=가속도, false=음압
  dbRef: 2.0e-5,               // dB 기준값 (음압: 20μPa)
);

final frequencies = result['frequencies'] as Float64List;
final magnitudes = result['magnitudes'] as Float64List; // dB
final amplitude = result['amplitude'] as Float64List;   // dB 변환 직전 선형 amplitude (RMS convention)
```

### Magnitude Spectrum

윈도우 함수 없이 전처리된 입력에 대해 순수 FFT를 수행합니다. 캘리브레이션 등 특수 용도에 사용됩니다.

```dart
final result = StftFftPlugin.computeMagnitudeSpectrum(
  preProcessedInput: calibrationData,
  samplingFrequency: 44100.0,
);
```

### STFT 스펙트로그램

시간-주파수 스펙트로그램을 계산하고, 색상 LUT를 적용하여 RGBA 픽셀 버퍼를 생성합니다.

```dart
final result = StftFftPlugin.computeStftSpectrogram(
  input: audioSamples,
  chunkSize: 1024,
  samplingFrequency: 44100.0,
  isAcc: false,
  dbRef: 2.0e-5,
  weightingType: 0,            // 0=없음, 1=A, 2=B, 3=C 가중치
  overlapPercent: 50.0,        // 윈도우 중첩률 (%) — 50일 때 hop=chunk/2
  rmsRangeMultiplier: 1.0,     // RMS 스무싱 범위 (1=스무싱 없음)
  targetWidth: 800,            // 출력 픽셀 너비
  targetHeight: 400,           // 출력 픽셀 높이
  startFreqHz: 20.0,
  endFreqHz: 20000.0,
  startTime: 0.0,
  endTime: 10.0,
  totalDuration: 10.0,
  colorLut: colorLookupTable,  // Int32List (256개 ARGB 값)
  autoCalcMinMax: true,
  minLevel: 0.0,
  maxLevel: 0.0,
);

final spectrogram = result['spectrogram'] as Float64List; // dB 스펙트로그램
final pixels = result['pixels'] as Uint8List;             // RGBA 픽셀 버퍼
final pixelWidth = result['pixelWidth'] as int;
final pixelHeight = result['pixelHeight'] as int;
final amplitude = result['amplitude'] as Float64List;     // raw grid(ampRows×ampCols) 선형 amplitude
final ampRows = result['ampRows'] as int;                 // = num_frames
final ampCols = result['ampCols'] as int;                 // = num_bins
```

### 필터링 및 픽셀 재생성

이미 계산된 스펙트로그램에서 주파수/시간 범위를 변경하거나 픽셀만 재생성할 때 사용합니다. 전체 STFT를 다시 계산하지 않으므로 효율적입니다.

```dart
// 스펙트로그램 필터링 + 픽셀 생성
final filtered = StftFftPlugin.filterAndGeneratePixels(
  spectrogram: spectrogram,
  specRows: specRows,
  specCols: specCols,
  samplingFrequency: 44100.0,
  chunkSize: 1024,
  totalDuration: 10.0,
  startFreqHz: 100.0,         // 변경된 주파수 범위
  endFreqHz: 8000.0,
  startTime: 2.0,             // 변경된 시간 범위
  endTime: 8.0,
  targetWidth: 800,
  targetHeight: 400,
  colorLut: colorLookupTable,
  minLevel: -80.0,
  maxLevel: 0.0,
);

// 색상 맵만 변경할 때 (필터링 불필요)
final pixels = StftFftPlugin.generatePixelBuffer(
  filteredData: filteredData,
  filteredRows: filteredRows,
  filteredCols: filteredCols,
  targetWidth: 800,
  targetHeight: 400,
  colorLut: newColorLut,
  minLevel: -80.0,
  maxLevel: 0.0,
);
```

## 아키텍처

```
┌─────────────────────────────────────────────┐
│  Dart (StftFftPlugin)                       │
│  └─ FFI 래퍼 · 메모리 관리 · 타입 변환       │
├─────────────────────────────────────────────┤
│  C Native (stft_fft_plugin)                 │
│  ├─ fft_engine.c      — FFT + Hanning 윈도우 │
│  ├─ welch_method.c    — Welch + 주파수 가중치 │
│  └─ stft_pipeline.c   — STFT + 리샘플링 + 픽셀│
├─────────────────────────────────────────────┤
│  KissFFT (kiss_fft.c, kiss_fftr.c)          │
│  └─ Double precision Real FFT               │
└─────────────────────────────────────────────┘
```

## 프로젝트 구조

```
flutter_stft_fft_plugin/
├── lib/
│   ├── stft_fft_plugin.dart               # 공개 API (StftFftPlugin 클래스)
│   └── src/
│       └── stft_fft_plugin_bindings.dart   # ffigen 자동생성 바인딩
├── src/
│   ├── stft_fft_plugin.h                  # C 헤더 (구조체 · 함수 정의)
│   ├── stft_fft_plugin.c                  # 버전 정보
│   ├── fft_engine.c                       # FFT 엔진
│   ├── welch_method.c                     # Welch's Method
│   ├── stft_pipeline.c                    # STFT 파이프라인
│   ├── signal_helpers.c                   # 주파수 가중치 · RMS 스무싱 공유 헬퍼
│   ├── signal_helpers.h                   # 공유 헬퍼 선언
│   ├── CMakeLists.txt                     # Android/데스크톱 빌드
│   └── kissfft/                           # KissFFT 라이브러리
├── ios/
│   └── flutter_stft_fft_plugin.podspec    # iOS 빌드 설정
├── android/
│   └── build.gradle                       # Android 빌드 설정
├── test/                                  # C 단위 테스트 (CMake + ctest)
├── docs/                                  # amplitude 규약 · 테스트 결과
└── ffigen.yaml                            # FFI 바인딩 생성 설정
```

## 개발

### 요구사항

- Flutter SDK ≥ 3.3.0
- Dart SDK ≥ 3.4.0
- Xcode (iOS 빌드)
- Android NDK (Android 빌드)

### FFI 바인딩 재생성

C 헤더(`src/stft_fft_plugin.h`)를 수정한 후 바인딩을 재생성합니다:

```bash
dart run ffigen
```

### C 테스트 실행

```bash
cd test
gcc -o test_fft_engine test_fft_engine.c ../src/fft_engine.c ../src/kissfft/kiss_fft.c ../src/kissfft/kiss_fftr.c -I../src -I../src/kissfft -lm -Dkiss_fft_scalar=double
./test_fft_engine
```

## 기술 스택

| 구분 | 기술 |
|------|------|
| 앱 프레임워크 | Flutter 3.3.0+ |
| 언어 바인딩 | Dart FFI |
| 바인딩 생성 | ffigen |
| 네이티브 코어 | C (C99) |
| FFT 라이브러리 | KissFFT (double precision) |
| iOS 빌드 | CocoaPods |
| Android 빌드 | CMake + Gradle |

## 라이선스

BSD 3-Clause License — Copyright (c) 2026, Eric.

KissFFT는 BSD 3-Clause License로 배포됩니다 (Copyright © 2003-2010 Mark Borgerding).

자세한 내용은 [LICENSE](LICENSE) 파일을 참고하세요.
