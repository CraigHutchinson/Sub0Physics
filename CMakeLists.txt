cmake_minimum_required(VERSION 3.21)

# Detect if Sub0Physics is being included inside another project
if(DEFINED PROJECT_NAME)
    set(IS_SUBPROJECT ON)
endif()

project(Sub0Physics VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)

option(SUB0PHYSICS_BUILD_TESTING    "Build unit-tests"                   ON)
option(SUB0PHYSICS_BUILD_BENCHMARKS "Build benchmarks"                   OFF)
option(SUB0PHYSICS_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(SUB0PHYSICS_ENABLE_ITT "Intel ITT task markers for VTune (fetches intel/ittapi)" OFF)

include(cmake/Dependencies.cmake)
# Benchmark, budget, profiler and vectorization tooling for the Sub0 family, from Sub0Math.
include("${Sub0Math_SOURCE_DIR}/cmake/Sub0MathTooling.cmake")

# ── Core library ─────────────────────────────────────────────────────────────

add_library(Sub0Physics STATIC
    src/bounded_periodic_tridiagonal.cpp
    src/circle_clearance.cpp
    src/convex_terrain.cpp
    src/corridor_builder.cpp
    src/damped_oscillator.cpp
    src/fidelity_selector.cpp
    src/height_extract.cpp
    src/height_page.cpp
    src/model_builder.cpp
    src/multibody_constraints.cpp
    src/multibody_dynamics.cpp
    src/multibody_elements.cpp
    src/multibody_kinematics.cpp
    src/multibody_ring.cpp
    src/polyline_profile.cpp
    src/prepared_profile.cpp
)

add_library(Sub0Physics::Sub0Physics ALIAS Sub0Physics)
# Determinism (docs/MULTIBODY.md section 8): no floating-point contraction anywhere in the library, and
# PUBLIC so consumers' translation units compile the inline headers (integrators, laws, heightAt) the same way.
sub0math_target_deterministic_fp(Sub0Physics PUBLIC)

target_compile_features(Sub0Physics PUBLIC cxx_std_23)
target_link_libraries(Sub0Physics PUBLIC Sub0Math::Sub0Math)

target_compile_options(Sub0Physics PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /Zc:preprocessor>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wdouble-promotion
        -fstrict-aliasing -Wstrict-aliasing=2>
    $<$<AND:$<BOOL:${SUB0PHYSICS_WARNINGS_AS_ERRORS}>,$<CXX_COMPILER_ID:MSVC>>:/WX>
    $<$<AND:$<BOOL:${SUB0PHYSICS_WARNINGS_AS_ERRORS}>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-Werror>
)

sub0math_target_perf_options(Sub0Physics)

# Profiler task markers (sub0physics/profiling/profile_task.hpp); compiled out unless enabled.
if(SUB0PHYSICS_ENABLE_ITT)
    CPMAddPackage(NAME ittapi GITHUB_REPOSITORY intel/ittapi
        GIT_TAG 69521ad9163f62b5e650f01c2e70f91595ec405e # v3.28.2
        EXCLUDE_FROM_ALL YES OPTIONS "ITT_API_INSTALL OFF")
    target_link_libraries(Sub0Physics PUBLIC ittnotify)
    target_compile_definitions(Sub0Physics PUBLIC SUB0PHYSICS_ITT=1)
endif()

target_include_directories(Sub0Physics
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(SUB0PHYSICS_BUILD_TESTING AND NOT IS_SUBPROJECT)
    enable_testing()
    add_subdirectory(tests)
endif()

# ── Install ───────────────────────────────────────────────────────────────────

# Installing needs Sub0Math and ittapi exported too; consumers use CPM instead until that is needed.
if(FALSE)

    include(CMakePackageConfigHelpers)

    set(SUB0PHYSICS_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Sub0Physics")

    configure_package_config_file(
        "${CMAKE_CURRENT_LIST_DIR}/cmake/Sub0PhysicsConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/Sub0PhysicsConfig.cmake"
        INSTALL_DESTINATION
            "${SUB0PHYSICS_CMAKE_CONFIG_DESTINATION}"
    )

    write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/Sub0PhysicsConfigVersion.cmake"
        COMPATIBILITY
            SameMinorVersion
    )

    install(
        TARGETS
            Sub0Physics
        EXPORT
            Sub0PhysicsTargets
        DESTINATION
            ${CMAKE_INSTALL_LIBDIR}
    )

    install(
        EXPORT
            Sub0PhysicsTargets
        NAMESPACE
            Sub0Physics::
        DESTINATION
            ${SUB0PHYSICS_CMAKE_CONFIG_DESTINATION}
    )

    install(
        FILES
            "${PROJECT_BINARY_DIR}/Sub0PhysicsConfigVersion.cmake"
            "${PROJECT_BINARY_DIR}/Sub0PhysicsConfig.cmake"
        DESTINATION
            "${SUB0PHYSICS_CMAKE_CONFIG_DESTINATION}"
    )

    install(
        DIRECTORY
            "include/"
        DESTINATION
            "${CMAKE_INSTALL_INCLUDEDIR}"
    )

endif()
