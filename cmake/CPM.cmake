# Verified bootstrap; keep this version consistent across sub0 consumers.
# 0.43+ reads CPM_<Package>_SOURCE from the environment, so a developer machine can point every
# package at a local checkout (CPM warns for each); CI sets none and builds the pins. Force a pin
# locally with -DCPM_<Package>_SOURCE= (defined and empty).
include_guard(GLOBAL)
set(CPM_DOWNLOAD_VERSION 0.43.2)
set(CPM_DOWNLOAD_SHA256 49a3bef91ceb65bb66d57255e12d1ffd22abc2f6408fa9fe4534c544a2f232aa)
set(_cpm_file "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
if(NOT EXISTS "${_cpm_file}")
    file(DOWNLOAD
        "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake"
        "${_cpm_file}" EXPECTED_HASH SHA256=${CPM_DOWNLOAD_SHA256}
        TLS_VERIFY ON STATUS _cpm_status)
    list(GET _cpm_status 0 _cpm_code)
    if(NOT _cpm_code EQUAL 0)
        file(REMOVE "${_cpm_file}")
        message(FATAL_ERROR "CPM download failed: ${_cpm_status}")
    endif()
endif()
file(SHA256 "${_cpm_file}" _cpm_hash)
if(NOT _cpm_hash STREQUAL CPM_DOWNLOAD_SHA256)
    message(FATAL_ERROR "CPM checksum mismatch: remove ${_cpm_file} and configure again.")
endif()
include("${_cpm_file}")
