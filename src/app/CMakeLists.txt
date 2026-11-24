# The MORPH driver: main.f90 with its input (config/), its mesh I/O and its
# observer of the run (io/). These are compiled into the executable only, so
# ORION and FiNeR never reach MORPHL.
#
# The reader of the MORPH boundary-condition file (io/Morph_BC_IO.f90) uses
# nothing but MORPHL. It is a library of its own, MORPHBC, which the driver and
# the unit programs that read such a file (test/) link: one target alone builds
# it and writes its module.
file(GLOB MORPH_APP_SRC CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/config/*.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/*.f90
)
list(REMOVE_ITEM MORPH_APP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/io/Morph_BC_IO.f90)

add_library(MORPHBC STATIC io/Morph_BC_IO.f90)
target_link_libraries(MORPHBC PUBLIC MORPHL)

add_executable(${solver} main.f90 ${MORPH_APP_SRC})
set_target_properties(${solver} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${solver} PRIVATE MORPHBC MORPHL ORION FiNeR::FiNeR)

if (USE_TECIO)
  target_link_libraries(${solver} PRIVATE -pthread)
endif()

if (USE_OPENMP)
  target_compile_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
endif()

set_target_properties(${solver} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
