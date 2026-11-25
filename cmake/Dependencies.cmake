include_guard(GLOBAL)
include("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")

# Pins, not HEADs: change deliberately and record why in the commit message.
# Sub0Math 67ede5a: the performance tooling moved out to Sub0Perf; sub0math_target_deterministic_fp is now Sub0Math's cmake/Sub0MathFp.cmake.
set(SUB0PHYSICS_SUB0MATH_REVISION 67ede5a2af8e17ade8ded1aaf008aad5ab43a3e9)
# Sub0Perf v0.4.0 (af59534, 2026-11-25): the family's performance tooling (ADR 0004 in CraigHutchinson/.github), fetched for its tooling only.
set(SUB0PHYSICS_SUB0PERF_REVISION af5953411488674d8f9b85e9ce8687510203c0ac)

CPMAddPackage(NAME Sub0Math GITHUB_REPOSITORY CraigHutchinson/Sub0Math
    GIT_TAG ${SUB0PHYSICS_SUB0MATH_REVISION} EXCLUDE_FROM_ALL YES
    OPTIONS "SUB0MATH_BUILD_TESTING OFF" "SUB0MATH_BUILD_BENCHMARKS OFF")

CPMAddPackage(NAME Sub0Perf GITHUB_REPOSITORY CraigHutchinson/Sub0Perf
    GIT_TAG ${SUB0PHYSICS_SUB0PERF_REVISION} DOWNLOAD_ONLY YES)
