# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(INTERVAL_CCDDownloadExternal)


# if(TIGHT_INCLUSION_WITH_GMP)
#   #GMP
#   find_package(GMPECCD)
#   IF(NOT ${GMP_FOUND})
#           MESSAGE(FATAL_ERROR "Cannot find GMP")
#   ENDIF()
# endif()


if(NOT TARGET Eigen3::Eigen)
  ccd_download_eigen()
  add_library(iccd_eigen INTERFACE)
  target_include_directories(iccd_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${INTERVAL_CCD_EXTERNAL}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET iccd_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS iccd_eigen)
endif()