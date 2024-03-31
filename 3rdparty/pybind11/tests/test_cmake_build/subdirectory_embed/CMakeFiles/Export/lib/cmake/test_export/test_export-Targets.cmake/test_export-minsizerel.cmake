#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "test_embed_lib" for configuration "MinSizeRel"
set_property(TARGET test_embed_lib APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(test_embed_lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/libtest_embed_lib.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS test_embed_lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_test_embed_lib "${_IMPORT_PREFIX}/bin/libtest_embed_lib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
