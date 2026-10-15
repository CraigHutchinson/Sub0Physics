# Performance tooling shared by Sub0Physics and its consumers.
#
#   sub0physics_target_perf_options(<target>)
#       Applies the profiling, vectorization-report and native-architecture options to a target.
#   sub0physics_add_benchmark(<target> <sources>...)
#       Builds a nanobench executable using tools/perf/bench_runner.hpp, plus these targets:
#         <target>_json             run it and write build/perf/<target>.json
#         <target>_check            compare with perf/baselines/<host>.json (fails on regression)
#         <target>_baseline         record perf/baselines/<host>.json from this host
#         <target>_vtune_hotspots   VTune hotspots      (SUB0PHYSICS_PROFILE_FILTER, _SECONDS)
#         <target>_vtune_uarch      VTune microarchitecture exploration
#         <target>_advisor_roofline Intel Advisor vectorization survey and roofline
#         <target>_perf_record      Linux perf sampling with call graphs
#
# Consumers include this file from the Sub0Physics source tree (CPM sets Sub0Physics_SOURCE_DIR).

include_guard(GLOBAL)

get_filename_component(SUB0PHYSICS_TOOLING_DIR "${CMAKE_CURRENT_LIST_DIR}/../tools/perf" ABSOLUTE)
get_filename_component(SUB0PHYSICS_BENCH_VENDOR_DIR "${CMAKE_CURRENT_LIST_DIR}/../tests/vendor" ABSOLUTE)

option(SUB0PHYSICS_PROFILE_FLAGS "Debug info and frame pointers in optimised builds, for profilers" OFF)
option(SUB0PHYSICS_VECTORIZATION_REPORT "Emit compiler vectorization remarks" OFF)
option(SUB0PHYSICS_NATIVE_ARCH "Tune for the build machine's CPU (local experiments only)" OFF)
set(SUB0PHYSICS_PROFILE_FILTER "" CACHE STRING "Benchmark name substring the profiling targets run")
set(SUB0PHYSICS_PROFILE_SECONDS "10" CACHE STRING "Soak seconds per benchmark under a profiler")

cmake_host_system_information(RESULT _sub0physics_host QUERY HOSTNAME)
string(TOLOWER "${_sub0physics_host}" _sub0physics_host)
set(SUB0PHYSICS_PERF_HOST "${_sub0physics_host}" CACHE STRING "Baseline name in perf/baselines/")

function(sub0physics_target_perf_options target)
    if(SUB0PHYSICS_PROFILE_FLAGS)
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/Zi /Oy->
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-g -fno-omit-frame-pointer>)
        target_link_options(${target} PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/DEBUG>)
    endif()
    if(SUB0PHYSICS_VECTORIZATION_REPORT)
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/Qvec-report:2>
            $<$<CXX_COMPILER_ID:GNU>:-fopt-info-vec-optimized -fopt-info-vec-missed>
            $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Rpass=loop-vectorize|slp-vectorizer
                -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize>)
    endif()
    if(SUB0PHYSICS_NATIVE_ARCH)
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-march=native>)
    endif()
endfunction()

function(sub0physics_add_benchmark target)
    add_executable(${target} ${ARGN})
    target_include_directories(${target} SYSTEM PRIVATE "${SUB0PHYSICS_TOOLING_DIR}" "${SUB0PHYSICS_BENCH_VENDOR_DIR}")
    target_compile_features(${target} PRIVATE cxx_std_23)
    sub0physics_target_perf_options(${target})

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(NOT Python3_Interpreter_FOUND)
        return()
    endif()
    set(results "${CMAKE_BINARY_DIR}/perf/${target}.json")
    set(baseline "${PROJECT_SOURCE_DIR}/perf/baselines/${SUB0PHYSICS_PERF_HOST}.json")
    set(check "${SUB0PHYSICS_TOOLING_DIR}/bench_check.py")
    set(profile "${SUB0PHYSICS_TOOLING_DIR}/profile.py")

    add_custom_target(${target}_json
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/perf"
        COMMAND $<TARGET_FILE:${target}> --json "${results}"
        DEPENDS ${target} USES_TERMINAL VERBATIM)
    add_custom_target(${target}_check
        COMMAND ${Python3_EXECUTABLE} "${check}" check --results "${results}" --baseline "${baseline}"
        DEPENDS ${target}_json USES_TERMINAL VERBATIM)
    add_custom_target(${target}_baseline
        COMMAND ${Python3_EXECUTABLE} "${check}" baseline --results "${results}" --baseline "${baseline}"
        DEPENDS ${target}_json USES_TERMINAL VERBATIM)

    foreach(analysis vtune-hotspots vtune-uarch advisor-roofline perf-record)
        string(REPLACE "-" "_" suffix "${analysis}")
        add_custom_target(${target}_${suffix}
            COMMAND ${Python3_EXECUTABLE} "${profile}" ${analysis} --exe $<TARGET_FILE:${target}>
                    --filter "${SUB0PHYSICS_PROFILE_FILTER}" --seconds ${SUB0PHYSICS_PROFILE_SECONDS}
                    --out "${CMAKE_BINARY_DIR}/profiles"
            DEPENDS ${target} USES_TERMINAL VERBATIM)
    endforeach()
endfunction()
