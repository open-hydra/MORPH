# MORPHL: every .f90 under src/lib. It links nothing: the library uses no
# module but iso_fortran_env and its own (test/run_morph_guards.sh, guard
# lib-isolated), so a solver adds it with add_subdirectory and links it, and
# MORPH's tests and driver build on it like any other consumer.
file(GLOB_RECURSE MORPH_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.f90
)

add_library(MORPHL ${MORPH_SRC})
add_library(MORPHL::MORPHL ALIAS MORPHL)

# The .mod files go to the project's module directory when it has one (MORPH's
# own build, or a parent that sets CMAKE_Fortran_MODULE_DIRECTORY), else to a
# directory of this target; either way a target that links MORPHL finds them.
if(CMAKE_Fortran_MODULE_DIRECTORY)
  set(MORPH_MODULE_DIR "${CMAKE_Fortran_MODULE_DIRECTORY}")
else()
  set(MORPH_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules")
endif()
set_target_properties(MORPHL PROPERTIES Fortran_MODULE_DIRECTORY "${MORPH_MODULE_DIR}")
target_include_directories(MORPHL INTERFACE $<BUILD_INTERFACE:${MORPH_MODULE_DIR}>)

IF(USE_OPENMP)
    target_compile_options(MORPHL PRIVATE ${OpenMP_Fortran_FLAGS})
ENDIF(USE_OPENMP)
