include_guard(GLOBAL)
include("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")

# Pins, not HEADs: change deliberately and record why in the commit message.
set(SUB0PHYSICS_SUB0MATH_REVISION e7afd8c35a20e22fd0bd3f2cada0a838deb7eecf)

CPMAddPackage(NAME Sub0Math GITHUB_REPOSITORY CraigHutchinson/Sub0Math
    GIT_TAG ${SUB0PHYSICS_SUB0MATH_REVISION} EXCLUDE_FROM_ALL YES
    OPTIONS "SUB0MATH_BUILD_TESTING OFF" "SUB0MATH_BUILD_BENCHMARKS OFF")
