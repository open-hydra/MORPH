# MORPH's test suite:  ctest --test-dir <build> [-L unit|guard|driver] [-j 2]
#
#   unit    the three unit programs MORPH came with (test/unit), built on
#           MORPHL, at 1 and 8 threads, and the two outputs compared
#   guard   the architecture guards of run_morph_guards.sh, one test each, the
#           build-configuration check, and MORPHL embedded in another project
#   driver  bin/MORPH on the three-block mesh of driver/multiblock: one run per
#           motion law at 1 and 8 threads, gated on the per-step discrete-GCL
#           residual and checked node by node (verify.py), the two outputs
#           compared, and four runs the driver must refuse
#
# An 8-thread test asks ctest for 8 processors (PROCESSORS), so `ctest -j N`
# does not oversubscribe a shared machine. In a serial build (USE_OPENMP OFF)
# the 8-thread tests and the comparisons are not registered.

set(MORPH_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MORPH_RUN_DIR  ${CMAKE_CURRENT_BINARY_DIR})

if(USE_OPENMP)
  set(MORPH_THREADS 1 8)
else()
  set(MORPH_THREADS 1)
  message(STATUS "MORPH tests: serial build, the 8-thread tests and the 1-vs-8 comparisons are not registered")
endif()

# ---------------------------------------------------------------------------
# unit
# ---------------------------------------------------------------------------
foreach(prog test_morph_block test_morph_gcl test_morph_nodeshift)
  add_executable(${prog} unit/${prog}.f90)
  set_target_properties(${prog} PROPERTIES LINKER_LANGUAGE Fortran
                        RUNTIME_OUTPUT_DIRECTORY ${MORPH_RUN_DIR}/unit)
  target_link_libraries(${prog} PRIVATE MORPHL)
  if(USE_OPENMP)
    target_compile_options(${prog} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${prog} PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()

  foreach(nt ${MORPH_THREADS})
    set(t unit_${prog}_omp${nt})
    add_test(NAME ${t}
             COMMAND bash ${MORPH_TEST_DIR}/unit/run_unit.sh $<TARGET_FILE:${prog}> ${nt}
                          ${MORPH_RUN_DIR}/unit/${prog}.omp${nt}.out)
    set_tests_properties(${t} PROPERTIES LABELS unit PROCESSORS ${nt} FIXTURES_SETUP ${t}
                         FAIL_REGULAR_EXPRESSION "FAIL")
  endforeach()

  if(USE_OPENMP)
    add_test(NAME unit_${prog}_threads_identical
             COMMAND bash ${MORPH_TEST_DIR}/same_files.sh
                          ${MORPH_RUN_DIR}/unit/${prog}.omp1.out ${MORPH_RUN_DIR}/unit/${prog}.omp8.out)
    set_tests_properties(unit_${prog}_threads_identical PROPERTIES LABELS unit
                         FIXTURES_REQUIRED "unit_${prog}_omp1;unit_${prog}_omp8")
  endif()
endforeach()

# ---------------------------------------------------------------------------
# guard
# ---------------------------------------------------------------------------
foreach(g no-fuss no-stop no-io lib-isolated)
  add_test(NAME guard_${g} COMMAND bash ${MORPH_TEST_DIR}/run_morph_guards.sh ${g})
  set_tests_properties(guard_${g} PROPERTIES LABELS guard)
endforeach()

# The -Wall guards compile with GNU flags, and the driver's against this
# tree's modules, which only the compiler that wrote them can read.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_test(NAME guard_wall-clean
           COMMAND bash ${MORPH_TEST_DIR}/run_morph_guards.sh wall-clean)
  add_test(NAME guard_wall-clean-app
           COMMAND bash ${MORPH_TEST_DIR}/run_morph_guards.sh --build ${CMAKE_BINARY_DIR} wall-clean-app)
  set_tests_properties(guard_wall-clean guard_wall-clean-app PROPERTIES LABELS guard
                       ENVIRONMENT "FC=${CMAKE_Fortran_COMPILER}")
else()
  message(STATUS "MORPH tests: the -Wall guards need gfortran; not registered with ${CMAKE_Fortran_COMPILER_ID}")
endif()

# The tree install.sh built and wrote the preset of; any other tree: skipped
add_test(NAME guard_build-config
         COMMAND bash ${MORPH_TEST_DIR}/check_build_config.sh --tree ${CMAKE_BINARY_DIR})
set_tests_properties(guard_build-config PROPERTIES LABELS guard SKIP_RETURN_CODE 77)

# What FUSS's switch to MORPH will do: add_subdirectory + link MORPHL, from a
# parent with a module directory and from one without (embed/run_embed.sh)
add_test(NAME guard_embed
         COMMAND bash ${MORPH_TEST_DIR}/embed/run_embed.sh ${CMAKE_SOURCE_DIR} ${MORPH_RUN_DIR}/embed
                      ${CMAKE_Fortran_COMPILER})
set_tests_properties(guard_embed PROPERTIES LABELS guard)

# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
find_package(Python3 3.6 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
  message(WARNING "MORPH tests: no python3 >= 3.6, so the driver tests (generate.py, verify.py) are not registered")
  return()
endif()

set(MORPH_CASE     ${MORPH_TEST_DIR}/driver/multiblock)
set(MORPH_RUN_CASE ${MORPH_TEST_DIR}/driver/run_case.sh)

foreach(law static translation prescribed nodeshift)
  foreach(nt ${MORPH_THREADS})
    set(t driver_${law}_omp${nt})
    add_test(NAME ${t}
             COMMAND bash ${MORPH_RUN_CASE} ${MORPH_CASE} input-${law}.ini
                          ${MORPH_RUN_DIR}/driver/${law}-omp${nt} ${nt} $<TARGET_FILE:${solver}>)
    set_tests_properties(${t} PROPERTIES LABELS driver PROCESSORS ${nt} FIXTURES_SETUP ${t}
                         ENVIRONMENT "PYTHON=${Python3_EXECUTABLE}")
  endforeach()

  if(USE_OPENMP)
    set(a ${MORPH_RUN_DIR}/driver/${law}-omp1/OUTPUT)
    set(b ${MORPH_RUN_DIR}/driver/${law}-omp8/OUTPUT)
    add_test(NAME driver_${law}_threads_identical
             COMMAND bash ${MORPH_TEST_DIR}/same_files.sh ${a}/mesh.p3d ${b}/mesh.p3d ${a}/gcl.dat ${b}/gcl.dat)
    set_tests_properties(driver_${law}_threads_identical PROPERTIES LABELS driver
                         FIXTURES_REQUIRED "driver_${law}_omp1;driver_${law}_omp8")
  endif()
endforeach()

# Runs the driver must refuse, one per exit status: 1 input, 2 a MORPH status,
# 3 the per-step GCL assertion (main.f90)
function(morph_refusal name ini code pattern)
  add_test(NAME driver_${name}
           COMMAND bash ${MORPH_RUN_CASE} ${MORPH_CASE} ${ini} ${MORPH_RUN_DIR}/driver/${name} 1
                        $<TARGET_FILE:${solver}> ${code} ${pattern})
  set_tests_properties(driver_${name} PROPERTIES LABELS driver ENVIRONMENT "PYTHON=${Python3_EXECUTABLE}")
endfunction()

morph_refusal(refuses_unknown_key   input-refuse-unknown.ini  1 "unknown option in input.ini: \\[MORPH-Motion\\] velocity")
morph_refusal(refuses_value_count   input-refuse-count.ini    1 "\\[MORPH-Motion\\] vel has 2 value\\(s\\), 3 expected")
morph_refusal(refuses_consumed_line input-refuse-consumed.ini 2 "failed in block 1, status [0-9]+ .*consumes a whole grid line")
morph_refusal(refuses_gcl_residual  input-refuse-gcl.ini      3 "step 1: discrete GCL residual .* exceeds gcl-tolerance")
