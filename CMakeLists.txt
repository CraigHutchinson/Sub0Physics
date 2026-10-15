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

# ── Core library ─────────────────────────────────────────────────────────────

add_library(Sub0Physics STATIC
    src/circle_clearance.cpp
    src/damped_oscillator.cpp
    src/fidelity_selector.cpp
    src/polyline_profile.cpp
)
add_library(Sub0Physics::Sub0Physics ALIAS Sub0Physics)

target_compile_features(Sub0Physics PUBLIC cxx_std_23)

target_compile_options(Sub0Physics PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /Zc:preprocessor>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wdouble-promotion>
    $<$<AND:$<BOOL:${SUB0PHYSICS_WARNINGS_AS_ERRORS}>,$<CXX_COMPILER_ID:MSVC>>:/WX>
    $<$<AND:$<BOOL:${SUB0PHYSICS_WARNINGS_AS_ERRORS}>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-Werror>
)

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

if(NOT IS_SUBPROJECT)

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
