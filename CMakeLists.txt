
cmake_minimum_required(VERSION 3.8)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)
if(STATIC_BUILD) # static build for linux system
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

project(RactIP CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE="ractip")
add_definitions(-DVERSION="1.1.1")
add_definitions(-DPACKAGE_VERSION="1.1.1")

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
elseif(ENABLE_SCIP)
  find_package(SCIP)
  if (SCIP_FOUND)
    MESSAGE(STATUS "Found SCIP version ${SCIP_VERSION}")
    add_definitions(-DWITH_SCIP)
    include_directories(${SCIP_INCLUDE_DIRS})
  endif()
elseif(ENABLE_HIGHS)
  find_package(HiGHS)
  find_package(ZLIB 1.2.3)
  if (HiGHS_FOUND)
    add_definitions(-DWITH_HIGHS)
    include_directories(${HiGHS_INCLUDE_DIR})
  endif()
else()
  find_package(GLPK)
  if (GLPK_FOUND)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
  endif()
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
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${SCIP_LIBRARIES} ${HiGHS_LIBRARY} ${ZLIB_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_features(ractip PUBLIC cxx_std_17)

install(TARGETS ractip RUNTIME DESTINATION bin)
