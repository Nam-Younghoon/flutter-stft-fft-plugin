cmake_minimum_required(VERSION 3.10)
project(stft_fft_plugin_tests LANGUAGES C)

enable_testing()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# 신규 — signal_helpers 단위 테스트 (ADR-002)
add_executable(signal_helpers_test_runner test_signal_helpers.c)
target_sources(signal_helpers_test_runner PRIVATE
  ${SRC_DIR}/fft_engine.c
  ${SRC_DIR}/welch_method.c
  ${SRC_DIR}/stft_pipeline.c
  ${SRC_DIR}/signal_helpers.c
  ${SRC_DIR}/kissfft/kiss_fft.c
  ${SRC_DIR}/kissfft/kiss_fftr.c
)
target_include_directories(signal_helpers_test_runner PRIVATE ${SRC_DIR} ${SRC_DIR}/kissfft)
target_compile_definitions(signal_helpers_test_runner PRIVATE kiss_fft_scalar=double)
target_link_libraries(signal_helpers_test_runner PRIVATE m)
add_test(NAME signal_helpers COMMAND signal_helpers_test_runner)

# 기존 회귀 가드 — Welch / STFT / FFT 출력이 변하지 않음
add_executable(welch_method_test_runner test_welch_method.c)
target_sources(welch_method_test_runner PRIVATE
  ${SRC_DIR}/fft_engine.c
  ${SRC_DIR}/welch_method.c
  ${SRC_DIR}/stft_pipeline.c
  ${SRC_DIR}/signal_helpers.c
  ${SRC_DIR}/kissfft/kiss_fft.c
  ${SRC_DIR}/kissfft/kiss_fftr.c
)
target_include_directories(welch_method_test_runner PRIVATE ${SRC_DIR} ${SRC_DIR}/kissfft)
target_compile_definitions(welch_method_test_runner PRIVATE kiss_fft_scalar=double)
target_link_libraries(welch_method_test_runner PRIVATE m)
add_test(NAME welch_method COMMAND welch_method_test_runner)

add_executable(stft_pipeline_test_runner test_stft_pipeline.c)
target_sources(stft_pipeline_test_runner PRIVATE
  ${SRC_DIR}/fft_engine.c
  ${SRC_DIR}/welch_method.c
  ${SRC_DIR}/stft_pipeline.c
  ${SRC_DIR}/signal_helpers.c
  ${SRC_DIR}/kissfft/kiss_fft.c
  ${SRC_DIR}/kissfft/kiss_fftr.c
)
target_include_directories(stft_pipeline_test_runner PRIVATE ${SRC_DIR} ${SRC_DIR}/kissfft)
target_compile_definitions(stft_pipeline_test_runner PRIVATE kiss_fft_scalar=double)
target_link_libraries(stft_pipeline_test_runner PRIVATE m)
add_test(NAME stft_pipeline COMMAND stft_pipeline_test_runner)

add_executable(fft_engine_test_runner test_fft_engine.c)
target_sources(fft_engine_test_runner PRIVATE
  ${SRC_DIR}/fft_engine.c
  ${SRC_DIR}/welch_method.c
  ${SRC_DIR}/stft_pipeline.c
  ${SRC_DIR}/signal_helpers.c
  ${SRC_DIR}/kissfft/kiss_fft.c
  ${SRC_DIR}/kissfft/kiss_fftr.c
)
target_include_directories(fft_engine_test_runner PRIVATE ${SRC_DIR} ${SRC_DIR}/kissfft)
target_compile_definitions(fft_engine_test_runner PRIVATE kiss_fft_scalar=double)
target_link_libraries(fft_engine_test_runner PRIVATE m)
add_test(NAME fft_engine COMMAND fft_engine_test_runner)
