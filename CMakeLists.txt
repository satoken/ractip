
cmake_minimum_required(VERSION 3.8)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)

project(RactIP CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE="ractip")
add_definitions(-DVERSION="1.0.2")
add_definitions(-DPACKAGE_VERSION="1.0.2")

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
include_directories(${VIENNARNA_INCLUDE_DIRS})
add_definitions(-DHAVE_VIENNA18)
add_definitions(-DHAVE_VIENNA20)

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
  endif()
else()
  if(STATIC_BUILD) # static build for linux system
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    find_library(GLPK_LIBRARY libglpk.a)
    find_library(ZLIB_LIBRARY libz.a)
    find_library(GMP_LIBRARY libgmp.a)
    find_library(DL_LIBRARY libdl.a)
    # find_library(LTDL_LIBRARY libltdl.a)
    # find_library(MPFR_LIBRARY libmpfr.a)
  else()
    find_library(GLPK_LIBRARY glpk)
    find_library(ZLIB_LIBRARY z)
    find_library(GMP_LIBRARY gmp)
    find_library(DL_LIBRARY dl)
    # find_library(LTDL_LIBRARY ltdl)
    # find_library(MPFR_LIBRARY mpfr)
  endif()
  add_definitions(-DWITH_GLPK)
endif()

add_executable(
  ractip
	src/ractip.cpp src/ip.cpp src/ip.h src/fa.cpp src/fa.h
	src/boltzmann_param.c src/boltzmann_param.h
	src/cmdline.c src/cmdline.h src/pf_duplex.c src/pf_duplex.h src/ushuffle.c src/ushuffle.h 
	src/contrafold/InferenceEngine.hpp src/contrafold/InferenceEngine.ipp 
	src/contrafold/DuplexEngine.hpp src/contrafold/DuplexEngine.ipp 
	src/contrafold/ParameterManager.hpp src/contrafold/ParameterManager.ipp 
	src/contrafold/SStruct.hpp src/contrafold/SStruct.cpp
	src/contrafold/Utilities.hpp src/contrafold/Utilities.ipp src/contrafold/Utilities.cpp
	src/contrafold/Defaults.ipp src/contrafold/LogSpace.hpp src/contrafold/Config.hpp src/contrafold/rand.h 
)
target_link_libraries(ractip 
  ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${CMAKE_DL_LIBS})
# target_compile_features(ractip PUBLIC cxx_std_17)

install(TARGETS ractip RUNTIME DESTINATION bin)
