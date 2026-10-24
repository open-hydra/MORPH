# The MORPH driver: main.f90 with its input (config/) and mesh I/O (io/). These
# are compiled into the executable only, so ORION and FiNeR never reach MORPHL.
file(GLOB MORPH_APP_SRC CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/config/*.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/io/*.f90
)

add_executable(${solver} main.f90 ${MORPH_APP_SRC})
set_target_properties(${solver} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${solver} PRIVATE MORPHL ORION FiNeR::FiNeR)

if (USE_TECIO)
  target_link_libraries(${solver} PRIVATE -pthread)
endif()

if (USE_OPENMP)
  target_compile_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(${solver} PRIVATE ${OpenMP_Fortran_FLAGS})
endif()

set_target_properties(${solver} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
