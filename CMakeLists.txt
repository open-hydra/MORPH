# CMake project file for 'MORPH'
#
# Two targets (plan 09 section 2):
#   MORPHL   the mesh-motion library, src/lib: kernels and the domain type; it
#            links nothing and uses no module but iso_fortran_env and its own
#   MORPH    the standalone driver, src/app: reads a mesh, applies a motion law,
#            writes the mesh; the only part that uses ORION (I/O) and FiNeR (INI)
# A solver that wants MORPH's geometry adds this directory and links MORPHL;
# the driver, ORION and FiNeR are then not configured at all.

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver MORPH)

project (${solver} LANGUAGES Fortran)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # MORPH has no MPI: it is rank-agnostic and operates on whatever block it is
  # handed (plan 09 section 9).
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
  # ORION's own default is ON, and with TecIO ORION builds TecIO inside ITS
  # source tree (lib/TecIO/<flavour>-install-<compiler>). Declared here first so
  # ORION's option() finds the cache entry and keeps OFF unless asked.
  OPTION(USE_TECIO "Link ORION against TecIO (binary Tecplot files)" OFF)

  # ------------------------------------------------------------------------
  # LTO-aware archiver (GNU). Copied from FUSS 2943354: the block is FUSS's
  # CMakeLists.txt:44-71 and this comment its :21-43, with MORPH's names.
  #
  # SetFortranFlags.cmake adds "-flto" to the RELEASE flags for GNU. The object
  # files then carry LTO bytecode instead of ordinary symbols, and plain
  # binutils ar/ranlib cannot read them ("plugin needed to handle lto object").
  # The result is a libMORPHL.a with an EMPTY symbol index and undefined
  # __morph_*_MOD_* references when linking bin/MORPH -- after the build has
  # already reported most of its progress. gcc-ar/gcc-ranlib/gcc-nm are thin
  # wrappers that load the LTO plugin and index the archive correctly.
  #
  # Two reasons this lives HERE and not in SetFortranFlags.cmake:
  #   1. That file returns early when the flags are already cached, so on an
  #      incremental re-configure (install.sh compile / cmake --preset) the
  #      block would be skipped and CMAKE_AR would silently revert.
  #   2. project() loads CMakeFortranCompiler.cmake, which does a plain
  #      "set(CMAKE_AR /usr/bin/ar)". A normal variable shadows the cache entry
  #      of the same name, so BOTH must be set -- cache for persistence, normal
  #      to win in this scope and to propagate into add_subdirectory() for the
  #      ORION / FiNeR sub-projects.
  # gcc-ar is a superset of ar, so this is applied for GNU regardless of build
  # type rather than being made conditional on the flag string.
  # ------------------------------------------------------------------------
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_AR AND CMAKE_Fortran_COMPILER_RANLIB)
      # version-matched wrappers CMake already detected (e.g. gcc-ar-12)
      set(LTO_AR     "${CMAKE_Fortran_COMPILER_AR}")
      set(LTO_RANLIB "${CMAKE_Fortran_COMPILER_RANLIB}")
    else()
      find_program(LTO_AR     NAMES gcc-ar     DOC "LTO-aware ar")
      find_program(LTO_RANLIB NAMES gcc-ranlib DOC "LTO-aware ranlib")
    endif()
    find_program(LTO_NM NAMES gcc-nm DOC "LTO-aware nm")

    if(LTO_AR AND LTO_RANLIB)
      set(CMAKE_AR     "${LTO_AR}"     CACHE FILEPATH "LTO-aware archiver" FORCE)
      set(CMAKE_RANLIB "${LTO_RANLIB}" CACHE FILEPATH "LTO-aware ranlib"   FORCE)
      set(CMAKE_AR     "${LTO_AR}")
      set(CMAKE_RANLIB "${LTO_RANLIB}")
      if(LTO_NM)
        set(CMAKE_NM "${LTO_NM}" CACHE FILEPATH "LTO-aware nm" FORCE)
        set(CMAKE_NM "${LTO_NM}")
      endif()
      message(STATUS "GNU toolchain: archiving with ${LTO_AR} / ${LTO_RANLIB}")
    else()
      message(WARNING
        "No LTO-aware ar/ranlib found (gcc-ar / gcc-ranlib). If the RELEASE "
        "build enables -flto, the static library will get an empty symbol "
        "index and linking bin/MORPH will fail with undefined references.")
    endif()
  endif()

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING (the same file as FUSS's).
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
  # Locate and set parallelization libraries.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # ORION and FiNeR: install.sh passes --include-orion= / --include-finer= as
  # ORION_PATH / FINER_PATH. Without them the hydra layout is tried, MORPH being
  # hydra/src/lib/MORPH and the libraries hydra/lib/ORION and
  # hydra/lib/third_party/FiNeR. A relative path is taken from this directory.
  foreach(_lib ORION FINER)
    if(_lib STREQUAL "ORION")
      set(_hydra "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/ORION")
    else()
      set(_hydra "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/third_party/FiNeR")
    endif()
    if(NOT DEFINED ${_lib}_PATH OR "${${_lib}_PATH}" STREQUAL "")
      if(EXISTS "${_hydra}/CMakeLists.txt")
        get_filename_component(_hydra "${_hydra}" ABSOLUTE)
        set(${_lib}_PATH "${_hydra}" CACHE PATH "${_lib} source tree")
        message(STATUS "${_lib}_PATH not given: using hydra's ${_hydra}")
      else()
        message(FATAL_ERROR "${_lib}_PATH is not set and ${_hydra} does not exist: "
                            "pass -D${_lib}_PATH=<path> (install.sh build --include-orion= / --include-finer=)")
      endif()
    endif()
    get_filename_component(${_lib}_PATH "${${_lib}_PATH}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    if(NOT EXISTS "${${_lib}_PATH}/CMakeLists.txt")
      message(FATAL_ERROR "${_lib}_PATH = ${${_lib}_PATH} holds no CMakeLists.txt")
    endif()
  endforeach()

  # The binary dir of a subproject belongs to THIS build tree (plan 09 section 11;
  # FUSS plan 10 D9: a fixed path under the source tree made every second tree
  # write its ORION and FiNeR objects into build/).
  if (NOT TARGET ORION)
    set ( orionD "${ORION_PATH}" "${CMAKE_CURRENT_BINARY_DIR}/lib/ORION")
    add_subdirectory(${orionD})
  endif()

  # FiNeR (v2.2.0+ tracks none of its five dependencies; install.sh fills src/third_party/ at pinned commits)
  if (NOT TARGET FiNeR)
    foreach(_dep PENF FACE FLAP BeFoR64 StringiFor)
      if (NOT EXISTS "${FINER_PATH}/src/third_party/${_dep}/CMakeLists.txt")
        message(FATAL_ERROR "FiNeR dependency ${_dep} is missing under ${FINER_PATH}/src/third_party/: run ./install.sh update --include-finer=${FINER_PATH} (or build), which clones the five dependencies at their pinned commits")
      endif()
    endforeach()
    set ( finerD "${FINER_PATH}" "${CMAKE_CURRENT_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

# Build application and tests (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)

  # Test suite: ctest --test-dir <build> (see test/CMakeLists.txt)
  enable_testing()
  add_subdirectory(test)
endif()
