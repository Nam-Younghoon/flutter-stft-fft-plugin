cmake_minimum_required(VERSION 3.10)
project(stft_fft_plugin_library VERSION 0.0.1 LANGUAGES C)

add_library(stft_fft_plugin SHARED
  "stft_fft_plugin.c"
  "fft_engine.c"
  "welch_method.c"
  "stft_pipeline.c"
  "kissfft/kiss_fft.c"
  "kissfft/kiss_fftr.c"
)

target_include_directories(stft_fft_plugin PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/kissfft"
)

target_compile_definitions(stft_fft_plugin PRIVATE
  DART_SHARED_LIB
  kiss_fft_scalar=double
)

set_target_properties(stft_fft_plugin PROPERTIES
  PUBLIC_HEADER stft_fft_plugin.h
  OUTPUT_NAME "stft_fft_plugin"
  C_STANDARD 99
)

if (ANDROID)
  target_link_options(stft_fft_plugin PRIVATE "-Wl,-z,max-page-size=16384")
  target_link_libraries(stft_fft_plugin m)
endif()
