set(VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor")

# ── Unit tests ────────────────────────────────────────────────────────────────

add_executable(Sub0Physics_Tests
    test_main.cpp
    test_channel_table.cpp
    test_circle_clearance.cpp
    test_damped_oscillator.cpp
    test_elements.cpp
    test_energy_ledger.cpp
    test_error_statistics.cpp
    test_fidelity_selector.cpp
    test_integrators.cpp
    test_polyline_profile.cpp
)

find_package(Threads REQUIRED)
target_link_libraries(Sub0Physics_Tests PRIVATE Sub0Physics::Sub0Physics Threads::Threads)
target_include_directories(Sub0Physics_Tests SYSTEM PRIVATE "${VENDOR_DIR}")
target_compile_options(Sub0Physics_Tests PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wshadow -Wconversion -Wsign-conversion>
)

add_test(NAME Sub0Physics_Tests COMMAND Sub0Physics_Tests)

# Replacement allocation operators count every heap allocation, so they live in their own executable.
add_executable(Sub0Physics_AllocationTests test_main.cpp test_allocations.cpp)
target_link_libraries(Sub0Physics_AllocationTests PRIVATE Sub0Physics::Sub0Physics Threads::Threads)
target_include_directories(Sub0Physics_AllocationTests SYSTEM PRIVATE "${VENDOR_DIR}")
add_test(NAME Sub0Physics_AllocationTests COMMAND Sub0Physics_AllocationTests)

# ── Header self-containment ───────────────────────────────────────────────────
# One translation unit per public header, each including only that header (twice,
# to exercise its include guard). A header that leans on another's transitive
# includes fails this build.

file(GLOB_RECURSE SUB0PHYSICS_PUBLIC_HEADERS CONFIGURE_DEPENDS
    RELATIVE "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/include/sub0physics/*.hpp")
set(SUB0PHYSICS_HEADER_CHECK_SOURCES "")
foreach(header IN LISTS SUB0PHYSICS_PUBLIC_HEADERS)
    string(REPLACE ".hpp" ".cpp" check_name "${header}")
    set(check_source "${CMAKE_CURRENT_BINARY_DIR}/header_check/${check_name}")
    file(CONFIGURE OUTPUT "${check_source}" CONTENT
        "#include \"${header}\"
#include \"${header}\"
")
    list(APPEND SUB0PHYSICS_HEADER_CHECK_SOURCES "${check_source}")
endforeach()

add_library(Sub0Physics_HeaderCheck OBJECT ${SUB0PHYSICS_HEADER_CHECK_SOURCES})
target_link_libraries(Sub0Physics_HeaderCheck PRIVATE Sub0Physics::Sub0Physics)

# ── Benchmarks (not run by ctest — run manually) ──────────────────────────────

if(SUB0PHYSICS_BUILD_BENCHMARKS)
    sub0math_add_benchmark(Sub0Physics_Bench bench_physics.cpp)
    target_link_libraries(Sub0Physics_Bench PRIVATE Sub0Physics::Sub0Physics)
endif()
