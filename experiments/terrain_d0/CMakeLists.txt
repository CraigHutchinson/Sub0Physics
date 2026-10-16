# Terrain representation decision experiment (D0). Not part of the library.
add_executable(terrain_d0 main.cpp)
target_compile_features(terrain_d0 PRIVATE cxx_std_23)
target_compile_options(terrain_d0 PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /utf-8>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra>)
sub0math_target_perf_options(terrain_d0)
