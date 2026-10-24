# A solver that embeds MORPH (test/embed/run_embed.sh): add_subdirectory and
# link MORPHL, nothing else. MORPH_DIR is the MORPH tree; WITH_MODULE_DIR says
# whether this parent sets a module directory of its own (FUSS does).
cmake_minimum_required(VERSION 3.23)
project(morph_embed LANGUAGES Fortran)

if(WITH_MODULE_DIR)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
endif()

add_subdirectory(${MORPH_DIR} morph)

add_executable(probe probe.f90)
target_link_libraries(probe PRIVATE MORPHL)
