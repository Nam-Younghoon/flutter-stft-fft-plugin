import 'dart:ffi';
import 'dart:io';
import 'dart:typed_data';

import 'package:ffi/ffi.dart';

import 'src/stft_fft_plugin_bindings.dart';

export 'src/stft_fft_plugin_bindings.dart';

const String _libName = 'stft_fft_plugin';

final DynamicLibrary _dylib = () {
  if (Platform.isMacOS || Platform.isIOS) {
    return DynamicLibrary.open('$_libName.framework/$_libName');
  }
  if (Platform.isAndroid || Platform.isLinux) {
    return DynamicLibrary.open('lib$_libName.so');
  }
  if (Platform.isWindows) {
    return DynamicLibrary.open('$_libName.dll');
  }
  throw UnsupportedError('Unknown platform: ${Platform.operatingSystem}');
}();

final StftFftPluginBindings nativeBindings = StftFftPluginBindings(_dylib);

class StftFftPlugin {
  StftFftPlugin._();

  static int version() => nativeBindings.stft_fft_plugin_version();

  static Map<String, Float64List> computeRealFft({
    required Float64List samples,
    required double samplingFrequency,
  }) {
    final inputPtr = malloc<Double>(samples.length);
    inputPtr.asTypedList(samples.length).setAll(0, samples);

    final result = nativeBindings.compute_real_fft(
      inputPtr,
      samples.length,
      samplingFrequency,
    );

    malloc.free(inputPtr);

    if (result == nullptr) {
      return {
        'frequencies': Float64List(0),
        'magnitudes': Float64List(0),
      };
    }

    final binCount = result.ref.bin_count;
    final frequencies = Float64List(binCount);
    final magnitudes = Float64List(binCount);

    final freqList = result.ref.frequencies.asTypedList(binCount);
    final magList = result.ref.magnitudes.asTypedList(binCount);
    frequencies.setAll(0, freqList);
    magnitudes.setAll(0, magList);

    nativeBindings.free_fft_result(result);

    return {
      'frequencies': frequencies,
      'magnitudes': magnitudes,
    };
  }

  static Map<String, dynamic> computeWelchSpectrum({
    required Float64List samples,
    required double samplingFrequency,
    required double frequencyResolution,
    required double overlapPercent,
    required int weightingType,
    required double rmsRangeMultiplier,
    required bool isAcc,
    required double dbRef,
  }) {
    final inputPtr = malloc<Double>(samples.length);
    inputPtr.asTypedList(samples.length).setAll(0, samples);

    final result = nativeBindings.compute_welch_spectrum(
      inputPtr,
      samples.length,
      samplingFrequency,
      frequencyResolution,
      overlapPercent,
      weightingType,
      rmsRangeMultiplier,
      isAcc ? 1 : 0,
      dbRef,
    );

    malloc.free(inputPtr);

    if (result == nullptr) {
      return {
        'frequencies': Float64List(0),
        'magnitudes': Float64List(0),
      };
    }

    final binCount = result.ref.bin_count;
    final frequencies = Float64List(binCount);
    final magnitudes = Float64List(binCount);

    final freqList = result.ref.frequencies.asTypedList(binCount);
    final magList = result.ref.magnitudes.asTypedList(binCount);
    frequencies.setAll(0, freqList);
    magnitudes.setAll(0, magList);

    nativeBindings.free_welch_result(result);

    return {
      'frequencies': frequencies,
      'magnitudes': magnitudes,
    };
  }

  static Map<String, Float64List> computeMagnitudeSpectrum({
    required Float64List preProcessedInput,
    required double samplingFrequency,
  }) {
    final inputPtr = malloc<Double>(preProcessedInput.length);
    inputPtr.asTypedList(preProcessedInput.length).setAll(0, preProcessedInput);

    final result = nativeBindings.compute_magnitude_spectrum(
      inputPtr,
      preProcessedInput.length,
      samplingFrequency,
    );

    malloc.free(inputPtr);

    if (result == nullptr) {
      return {
        'frequencies': Float64List(0),
        'magnitudes': Float64List(0),
      };
    }

    final binCount = result.ref.bin_count;
    final frequencies = Float64List(binCount);
    final magnitudes = Float64List(binCount);

    final freqList = result.ref.frequencies.asTypedList(binCount);
    final magList = result.ref.magnitudes.asTypedList(binCount);
    frequencies.setAll(0, freqList);
    magnitudes.setAll(0, magList);

    nativeBindings.free_fft_result(result);

    return {
      'frequencies': frequencies,
      'magnitudes': magnitudes,
    };
  }

  static Map<String, dynamic> computeStftSpectrogram({
    required Float64List input,
    required int chunkSize,
    required double samplingFrequency,
    required bool isAcc,
    required double dbRef,
    required int targetWidth,
    required int targetHeight,
    required double startFreqHz,
    required double endFreqHz,
    required double startTime,
    required double endTime,
    required double totalDuration,
    required Int32List colorLut,
    required bool autoCalcMinMax,
    required double minLevel,
    required double maxLevel,
  }) {
    final inputPtr = malloc<Double>(input.length);
    inputPtr.asTypedList(input.length).setAll(0, input);

    final lutPtr = malloc<Int32>(colorLut.length);
    lutPtr.asTypedList(colorLut.length).setAll(0, colorLut);

    final result = nativeBindings.compute_stft_spectrogram(
      inputPtr,
      input.length,
      chunkSize,
      samplingFrequency,
      isAcc ? 1 : 0,
      dbRef,
      targetWidth,
      targetHeight,
      startFreqHz,
      endFreqHz,
      startTime,
      endTime,
      totalDuration,
      lutPtr,
      colorLut.length,
      autoCalcMinMax ? 1 : 0,
      minLevel,
      maxLevel,
    );

    malloc.free(inputPtr);
    malloc.free(lutPtr);

    if (result == nullptr) {
      return {
        'spectrogram': Float64List(0),
        'specRows': 0,
        'specCols': 0,
        'pixels': Uint8List(0),
        'pixelWidth': 0,
        'pixelHeight': 0,
        'minLevel': 0.0,
        'maxLevel': 0.0,
      };
    }

    final specRows = result.ref.spec_rows;
    final specCols = result.ref.spec_cols;
    final specTotal = specRows * specCols;
    final spectrogram = Float64List(specTotal);
    spectrogram.setAll(0, result.ref.spectrogram.asTypedList(specTotal));

    final pixelWidth = result.ref.pixel_width;
    final pixelHeight = result.ref.pixel_height;
    final pixelTotal = pixelWidth * pixelHeight * 4;
    final pixels = Uint8List(pixelTotal);
    pixels.setAll(0, result.ref.pixels.asTypedList(pixelTotal));

    final minVal = result.ref.min_level;
    final maxVal = result.ref.max_level;

    nativeBindings.free_stft_result(result);

    return {
      'spectrogram': spectrogram,
      'specRows': specRows,
      'specCols': specCols,
      'pixels': pixels,
      'pixelWidth': pixelWidth,
      'pixelHeight': pixelHeight,
      'minLevel': minVal,
      'maxLevel': maxVal,
    };
  }

  static Map<String, dynamic> filterAndGeneratePixels({
    required Float64List spectrogram,
    required int specRows,
    required int specCols,
    required double samplingFrequency,
    required int chunkSize,
    required double totalDuration,
    required double startFreqHz,
    required double endFreqHz,
    required double startTime,
    required double endTime,
    required int targetWidth,
    required int targetHeight,
    required Int32List colorLut,
    required double minLevel,
    required double maxLevel,
  }) {
    final specPtr = malloc<Double>(spectrogram.length);
    specPtr.asTypedList(spectrogram.length).setAll(0, spectrogram);

    final lutPtr = malloc<Int32>(colorLut.length);
    lutPtr.asTypedList(colorLut.length).setAll(0, colorLut);

    final result = nativeBindings.filter_and_generate_pixels(
      specPtr,
      specRows,
      specCols,
      samplingFrequency,
      chunkSize,
      totalDuration,
      startFreqHz,
      endFreqHz,
      startTime,
      endTime,
      targetWidth,
      targetHeight,
      lutPtr,
      colorLut.length,
      minLevel,
      maxLevel,
    );

    malloc.free(specPtr);
    malloc.free(lutPtr);

    if (result == nullptr) {
      return {
        'filteredData': Float64List(0),
        'filteredRows': 0,
        'filteredCols': 0,
        'pixels': Uint8List(0),
        'pixelWidth': 0,
        'pixelHeight': 0,
      };
    }

    final filteredRows = result.ref.filtered_rows;
    final filteredCols = result.ref.filtered_cols;
    final filteredTotal = filteredRows * filteredCols;
    final filteredData = Float64List(filteredTotal);
    filteredData.setAll(0, result.ref.filtered_data.asTypedList(filteredTotal));

    final pixelWidth = result.ref.pixel_width;
    final pixelHeight = result.ref.pixel_height;
    final pixelTotal = pixelWidth * pixelHeight * 4;
    final pixels = Uint8List(pixelTotal);
    pixels.setAll(0, result.ref.pixels.asTypedList(pixelTotal));

    nativeBindings.free_filtered_pixel_result(result);

    return {
      'filteredData': filteredData,
      'filteredRows': filteredRows,
      'filteredCols': filteredCols,
      'pixels': pixels,
      'pixelWidth': pixelWidth,
      'pixelHeight': pixelHeight,
    };
  }

  static Map<String, dynamic> generatePixelBuffer({
    required Float64List filteredData,
    required int filteredRows,
    required int filteredCols,
    required int targetWidth,
    required int targetHeight,
    required Int32List colorLut,
    required double minLevel,
    required double maxLevel,
  }) {
    final dataPtr = malloc<Double>(filteredData.length);
    dataPtr.asTypedList(filteredData.length).setAll(0, filteredData);

    final lutPtr = malloc<Int32>(colorLut.length);
    lutPtr.asTypedList(colorLut.length).setAll(0, colorLut);

    final result = nativeBindings.generate_pixel_buffer(
      dataPtr,
      filteredRows,
      filteredCols,
      targetWidth,
      targetHeight,
      lutPtr,
      colorLut.length,
      minLevel,
      maxLevel,
    );

    malloc.free(dataPtr);
    malloc.free(lutPtr);

    if (result == nullptr) {
      return {
        'pixels': Uint8List(0),
        'pixelWidth': 0,
        'pixelHeight': 0,
      };
    }

    final pixelWidth = result.ref.pixel_width;
    final pixelHeight = result.ref.pixel_height;
    final pixelTotal = pixelWidth * pixelHeight * 4;
    final pixels = Uint8List(pixelTotal);
    pixels.setAll(0, result.ref.pixels.asTypedList(pixelTotal));

    nativeBindings.free_pixel_buffer_result(result);

    return {
      'pixels': pixels,
      'pixelWidth': pixelWidth,
      'pixelHeight': pixelHeight,
    };
  }
}
