Pod::Spec.new do |s|
  s.name             = 'flutter_stft_fft_plugin'
  s.version          = '0.0.1'
  s.summary          = 'KissFFT-based native FFT/STFT plugin for Flutter.'
  s.homepage         = 'https://github.com/eric/flutter_stft_fft_plugin'
  s.license          = { :file => '../LICENSE' }
  s.author           = { 'Eric' => 'eric@example.com' }
  s.source           = { :path => '.' }
  s.source_files     = 'Classes/**/*'
  s.dependency 'Flutter'
  s.platform         = :ios, '12.0'

  s.pod_target_xcconfig = {
    'DEFINES_MODULE' => 'YES',
    'EXCLUDED_ARCHS[sdk=iphonesimulator*]' => 'i386',
    'GCC_PREPROCESSOR_DEFINITIONS' => 'kiss_fft_scalar=double DART_SHARED_LIB=1',
  }
  s.swift_version = '5.0'
end
