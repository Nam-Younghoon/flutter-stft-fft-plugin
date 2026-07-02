# flutter_stft_fft_plugin

[한국어](README.md) | **English**

A native FFT/STFT Flutter plugin built on KissFFT. The high-performance signal processing routines are written in C and invoked directly through Dart FFI.

## Features

- **Double-precision FFT** — Accurate frequency analysis backed by the KissFFT library
- **Welch's Method** — Stable power spectrum estimation via overlapped window averaging
- **STFT spectrogram** — Time–frequency 2D analysis with real-time pixel rendering
- **A/B/C frequency weighting** — Compliant with international acoustic measurement standards
- **Direct Dart FFI calls** — Native-level performance without platform-channel overhead

## Supported platforms

| Platform | Minimum version | Build system |
|----------|-----------------|--------------|
| iOS | 12.0 | CocoaPods |
| Android | API 24 | CMake |

## Reference docs

- [STFT / FFT amplitude · dB convention](docs/STFT_FFT_AMPLITUDE_CONVENTION.md) — amplitude definitions, the dB conversion formula, the `isAcc` unit assumption, and where ENBW correction is applied; everything you need to compare result dB values against the ISO 1683 standard.
- [v0.0.2 test results](docs/test_results_v0.0.2.md) — regression/new test pass log after the Welch RMS-convention fix.

## Installation

Add the dependency to your `pubspec.yaml`:

```yaml
dependencies:
  flutter_stft_fft_plugin:
    git:
      url: https://github.com/Nam-Younghoon/flutter-stft-fft-plugin.git
      ref: main
```

## API

All methods are exposed as static members of the `StftFftPlugin` class.

### Real FFT

Transforms a time-domain signal into the frequency domain. A Hanning window is applied automatically.

```dart
import 'package:flutter_stft_fft_plugin/stft_fft_plugin.dart';

final result = StftFftPlugin.computeRealFft(
  samples: audioSamples,        // Float64List
  samplingFrequency: 44100.0,
);

final frequencies = result['frequencies']!; // Float64List — frequency axis (Hz)
final magnitudes = result['magnitudes']!;   // Float64List — magnitude
```

### Welch's Method

Estimates a noise-reduced power spectrum using overlapped window averaging. Results are returned on a dB scale.

```dart
final result = StftFftPlugin.computeWelchSpectrum(
  samples: audioSamples,
  samplingFrequency: 44100.0,
  frequencyResolution: 1.0,    // Frequency resolution in Hz
  overlapPercent: 50.0,        // Window overlap (%)
  weightingType: 1,            // 0=none, 1=A, 2=B, 3=C weighting
  rmsRangeMultiplier: 1.0,     // RMS smoothing range
  isAcc: false,                // true=acceleration, false=sound pressure
  dbRef: 2.0e-5,               // dB reference (sound pressure: 20μPa)
);

final frequencies = result['frequencies'] as Float64List;
final magnitudes = result['magnitudes'] as Float64List; // dB
final amplitude = result['amplitude'] as Float64List;   // linear amplitude before dB conversion (RMS convention)
```

### Magnitude Spectrum

Runs a plain FFT over an already pre-processed input without applying a window function. Useful for specialised cases such as calibration.

```dart
final result = StftFftPlugin.computeMagnitudeSpectrum(
  preProcessedInput: calibrationData,
  samplingFrequency: 44100.0,
);
```

### STFT Spectrogram

Computes a time–frequency spectrogram and produces an RGBA pixel buffer by applying a colour LUT.

```dart
final result = StftFftPlugin.computeStftSpectrogram(
  input: audioSamples,
  chunkSize: 1024,
  samplingFrequency: 44100.0,
  isAcc: false,
  dbRef: 2.0e-5,
  weightingType: 0,            // 0=none, 1=A, 2=B, 3=C weighting
  overlapPercent: 50.0,        // Window overlap (%) — hop=chunk/2 at 50
  rmsRangeMultiplier: 1.0,     // RMS smoothing range (1=no smoothing)
  targetWidth: 800,            // Output pixel width
  targetHeight: 400,           // Output pixel height
  startFreqHz: 20.0,
  endFreqHz: 20000.0,
  startTime: 0.0,
  endTime: 10.0,
  totalDuration: 10.0,
  colorLut: colorLookupTable,  // Int32List (256 ARGB entries)
  autoCalcMinMax: true,
  minLevel: 0.0,
  maxLevel: 0.0,
);

final spectrogram = result['spectrogram'] as Float64List; // dB spectrogram
final pixels = result['pixels'] as Uint8List;             // RGBA pixel buffer
final pixelWidth = result['pixelWidth'] as int;
final pixelHeight = result['pixelHeight'] as int;
final amplitude = result['amplitude'] as Float64List;     // linear amplitude on the raw grid (ampRows×ampCols)
final ampRows = result['ampRows'] as int;                 // = num_frames
final ampCols = result['ampCols'] as int;                 // = num_bins
```

### Filtering and pixel regeneration

Use these helpers to adjust the frequency/time range or regenerate pixels from a spectrogram that has already been computed. They avoid recomputing the full STFT, making them efficient for interactive views.

```dart
// Filter the spectrogram and produce new pixels
final filtered = StftFftPlugin.filterAndGeneratePixels(
  spectrogram: spectrogram,
  specRows: specRows,
  specCols: specCols,
  samplingFrequency: 44100.0,
  chunkSize: 1024,
  totalDuration: 10.0,
  startFreqHz: 100.0,         // Updated frequency range
  endFreqHz: 8000.0,
  startTime: 2.0,             // Updated time range
  endTime: 8.0,
  targetWidth: 800,
  targetHeight: 400,
  colorLut: colorLookupTable,
  minLevel: -80.0,
  maxLevel: 0.0,
);

// Only update the colour map (no filtering required)
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

## Architecture

```
┌─────────────────────────────────────────────┐
│  Dart (StftFftPlugin)                       │
│  └─ FFI wrapper · memory management · types │
├─────────────────────────────────────────────┤
│  C native (stft_fft_plugin)                 │
│  ├─ fft_engine.c      — FFT + Hanning window │
│  ├─ welch_method.c    — Welch + freq weights │
│  └─ stft_pipeline.c   — STFT + resample + px │
├─────────────────────────────────────────────┤
│  KissFFT (kiss_fft.c, kiss_fftr.c)          │
│  └─ Double precision real FFT               │
└─────────────────────────────────────────────┘
```

## Project layout

```
flutter_stft_fft_plugin/
├── lib/
│   ├── stft_fft_plugin.dart               # Public API (StftFftPlugin class)
│   └── src/
│       └── stft_fft_plugin_bindings.dart   # ffigen-generated bindings
├── src/
│   ├── stft_fft_plugin.h                  # C header (structs · function decls)
│   ├── stft_fft_plugin.c                  # Version info
│   ├── fft_engine.c                       # FFT engine
│   ├── welch_method.c                     # Welch's Method
│   ├── stft_pipeline.c                    # STFT pipeline
│   ├── signal_helpers.c                   # Shared frequency-weighting · RMS-smoothing helpers
│   ├── signal_helpers.h                   # Shared helper declarations
│   ├── CMakeLists.txt                     # Android / desktop build
│   └── kissfft/                           # KissFFT library
├── ios/
│   └── flutter_stft_fft_plugin.podspec    # iOS build configuration
├── android/
│   └── build.gradle                       # Android build configuration
├── test/                                  # C unit tests (CMake + ctest)
├── docs/                                  # amplitude convention · test results
└── ffigen.yaml                            # FFI binding generator config
```

## Development

### Requirements

- Flutter SDK ≥ 3.3.0
- Dart SDK ≥ 3.4.0
- Xcode (iOS build)
- Android NDK (Android build)

### Regenerating the FFI bindings

After editing the C header (`src/stft_fft_plugin.h`), regenerate the bindings:

```bash
dart run ffigen --config ffigen.yaml
```

### Running the C tests

```bash
cd test
cc -o test_fft_engine test_fft_engine.c \
   ../src/fft_engine.c ../src/stft_fft_plugin.c \
   ../src/kissfft/kiss_fft.c ../src/kissfft/kiss_fftr.c \
   -I../src -I../src/kissfft -lm -Dkiss_fft_scalar=double
./test_fft_engine
```

## Technology stack

| Area | Technology |
|------|------------|
| App framework | Flutter 3.3.0+ |
| Language bindings | Dart FFI |
| Binding generator | ffigen |
| Native core | C (C99) |
| FFT library | KissFFT (double precision) |
| iOS build | CocoaPods |
| Android build | CMake + Gradle |

## License

BSD 3-Clause License — Copyright (c) 2026, Eric.

KissFFT is distributed under the BSD 3-Clause License (Copyright © 2003-2010 Mark Borgerding).

See the [LICENSE](LICENSE) file for details.
