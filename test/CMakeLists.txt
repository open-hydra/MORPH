# MORPH's test suite:  ctest --test-dir <build> [-L unit|guard|driver] [-j 2]
#
#   unit    the unit programs (test/unit): the three MORPH came with, the
#           smoothed interior, the event list, the given law, the one-call run,
#           the surface, the reader of the MORPH boundary file, the
#           level-set kernels and the level-set grid of blocks with its
#           exchanges (its simulated ranks are threads, so it asks for 4
#           processors and is stopped after 300 s if a rank is left
#           waiting); built on MORPHL (the surface and the reader
#           tests also on the driver's reader, MORPHBC), at 1 and 8 threads,
#           and the two outputs compared
#   guard   the architecture guards of run_morph_guards.sh, one test each, the
#           build-configuration check, and MORPHL embedded in another project
#   driver  bin/MORPH on the three-block mesh of driver/multiblock: one run per
#           motion law and per nodeshift variant (interior line or smooth, a
#           uniform or a bumped face displacement) at 1 and 8 threads, gated on
#           the per-step discrete-GCL residual and checked node by node
#           (verify.py), the two outputs compared, and seven runs the driver
#           must refuse
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
foreach(prog test_morph_block test_morph_gcl test_morph_nodeshift test_morph_smooth test_morph_events
             test_morph_given test_morph_run test_morph_surface test_morph_bc test_morph_ls
             test_morph_ls_blocks)
  add_executable(${prog} unit/${prog}.f90)
  set_target_properties(${prog} PROPERTIES LINKER_LANGUAGE Fortran
                        RUNTIME_OUTPUT_DIRECTORY ${MORPH_RUN_DIR}/unit)
  target_link_libraries(${prog} PRIVATE MORPHL)
  # the two that read a MORPH boundary file use the driver's reader
  if(prog MATCHES "^test_morph_(surface|bc)$")
    target_link_libraries(${prog} PRIVATE MORPHBC)
  endif()
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
    if(prog STREQUAL "test_morph_ls_blocks")
      if(nt LESS 4)
        set_tests_properties(${t} PROPERTIES PROCESSORS 4)
      endif()
      set_tests_properties(${t} PROPERTIES TIMEOUT 300)
    endif()
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

# What FUSS does: add_subdirectory + link MORPHL, from a parent with a module
# directory and from one without (embed/run_embed.sh)
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

# nodeshift-smooth: the uniform recession with the smoothed interior (it must
# reproduce the line law); nodeshift-bump: a face displacement that varies
# across the face, smoothed (the equilibrium checks); nodeshift-bump-line: the
# same displacement under the line law (the replica)
#
# Besides verify.py's checks, run.log must show (run_case.sh, exit 0 and the
# patterns): the interior of a nodeshift run and, when it is smoothed, the
# smoothing's settings and its report under the step lines (src/app/main.f90,
# src/app/io/Morph_Output.f90)
set(MORPH_LINE_INTERIOR_LINE   "^ -  Interior +--> line$")
set(MORPH_LINE_INTERIOR_SMOOTH "^ -  Interior +--> smooth$")
set(MORPH_LINE_SMOOTHING       "^ -  Smoothing +--> tolerance +[0-9.]+E[-+][0-9]+, maxit [0-9]+, omega ")
set(MORPH_LINE_SMOOTH_REPORT   "^ +smoothing: sweeps [0-9]+ \\(faces\\) \\+ [0-9]+ \\(interior\\), last update ")
set(MORPH_SHOWS_nodeshift           ${MORPH_LINE_INTERIOR_LINE})
set(MORPH_SHOWS_nodeshift-bump-line ${MORPH_LINE_INTERIOR_LINE})
set(MORPH_SHOWS_nodeshift-smooth    ${MORPH_LINE_INTERIOR_SMOOTH} ${MORPH_LINE_SMOOTHING} ${MORPH_LINE_SMOOTH_REPORT})
set(MORPH_SHOWS_nodeshift-bump      ${MORPH_LINE_INTERIOR_SMOOTH} ${MORPH_LINE_SMOOTHING} ${MORPH_LINE_SMOOTH_REPORT})
foreach(law static translation prescribed nodeshift nodeshift-smooth nodeshift-bump nodeshift-bump-line)
  set(shows)
  if(DEFINED MORPH_SHOWS_${law})
    set(shows 0 ${MORPH_SHOWS_${law}})
  endif()
  foreach(nt ${MORPH_THREADS})
    set(t driver_${law}_omp${nt})
    add_test(NAME ${t}
             COMMAND bash ${MORPH_RUN_CASE} ${MORPH_CASE} input-${law}.ini
                          ${MORPH_RUN_DIR}/driver/${law}-omp${nt} ${nt} $<TARGET_FILE:${solver}> ${shows})
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

# Runs the driver must refuse, by exit status: 1 input (the motion parameters
# included, which Morph_Run refuses), 2 a MORPH status in a step (not in the
# setup: the pattern says which), 3 the per-step GCL tolerance of Morph_Run
function(morph_refusal name ini code pattern)
  add_test(NAME driver_${name}
           COMMAND bash ${MORPH_RUN_CASE} ${MORPH_CASE} ${ini} ${MORPH_RUN_DIR}/driver/${name} 1
                        $<TARGET_FILE:${solver}> ${code} ${pattern})
  set_tests_properties(driver_${name} PROPERTIES LABELS driver ENVIRONMENT "PYTHON=${Python3_EXECUTABLE}")
endfunction()

morph_refusal(refuses_unknown_key   input-refuse-unknown.ini  1 "unknown option in input.ini: \\[MORPH-Motion\\] velocity")
morph_refusal(refuses_value_count   input-refuse-count.ini    1 "\\[MORPH-Motion\\] vel has 2 value\\(s\\), 3 expected")
morph_refusal(refuses_shift_block   input-refuse-block.ini    1 "\\[MORPH-Motion\\] run: nodeshift block 4 is not a block of the domain")
morph_refusal(refuses_consumed_line input-refuse-consumed.ini 2 "MORPH step failed in block 1, status [0-9]+ .*consumes a whole grid line")
morph_refusal(refuses_smooth_omega  input-refuse-smooth-omega.ini 1 "\\[MORPH-Motion\\] run: smooth omega +2\\.000E\\+00 must be below 2")
morph_refusal(refuses_smooth_maxit  input-refuse-smooth-maxit.ini 2 "MORPH step failed in block 1, status 9 .*smoothing: not converged in 1 sweeps.*tolerance +1\\.00E-10")
morph_refusal(refuses_gcl_residual  input-refuse-gcl.ini      3 "step 1: discrete GCL residual .* exceeds gcl-tolerance")
