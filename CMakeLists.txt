cmake_minimum_required(VERSION 3.14)

# set the project name
set(CMAKE_PROJECT_NAME "TOMOATT")
project(${CMAKE_PROJECT_NAME} C CXX)

# Default to C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

#
# for production
#
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -O3 -funroll-loops -ffast-math -march=native -ftree-vectorize -fopenmp")
#
# for debugging
#
#set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -g -funroll-loops") # debug
#
# for vtune profiling
# with vectorization report
#set(gcc_coverage_compile_flags "-lm -g -O3") # debug
#set(GCC_COVERAGE_COMPILE_FLAGS "-g") # debug
#
# for FUGAKU (use clang mode)
#
#set(GCC_COVERAGE_COMPILE_FLAGS "-Nclang -Ofast -std=c++17")
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

# find installed MPI
message(STATUS "Running CMAKE FindMPI.cmake...")
find_package(MPI)

message(STATUS "MPI_FOUND: ${MPI_FOUND}")
message(STATUS "MPI_VERSION: ${MPI_VERSION}")

# find openmp
find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS "OpenMP found")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  add_definitions(-DUSE_OMP)
endif()

# find HDF5 parallel TODO: check parallel io is enable or not
message(STATUS "Running CMAKE FindHDF5.cmake...")
# set parallel HDF5 default
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5)
if(HDF5_FOUND)
  message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
  add_definitions(-DUSE_HDF5)
endif()

# find cuda package if USE_CUDA is defined
if(USE_CUDA)
  message(STATUS "Running CMAKE FindCUDA.cmake...")
  enable_language(CUDA)
  find_package(CUDA)
else()
  message(STATUS "TomoATT is compiled without cuda, because -DUSE_CUDA=True is not defined")
endif()

if(CUDA_FOUND) # TODO : add HIP here in the future
  message(STATUS "CUDA_FOUND: ${CUDA_FOUND}")
  add_definitions(-DUSE_CUDA)
  set(CUDA_LIBRARY_NAME "TOMOATT_CUDA")

  # list of source and header files for cuda
  file(GLOB SOURCES_CUDA "cuda/*.cu")
  file(GLOB HEADERS_CUDA "cuda/*.cuh")

  # cuda flag
  #
  # for production
  set(CMAKE_CUDA_FLAGS "-fPIC -O3 -use_fast_math -extra-device-vectorization  -gencode arch=compute_61,code=sm_61")
  #
  # for debugging
  #set(CMAKE_CUDA_FLAGS "-fPIC -lineinfo -g -O0 -gencode arch=compute_61,code=sm_61")
  set(CMAKE_CUDA_STANDARD "11")
  message(STATUS, "TomoATT will be compiled with cuda.")
else()
  message(STATUS "TomoATT will be compiled without cuda, because cuda is not found or -DUSE_CUDA=True was not specified.")
endif()


# submodules
# yaml parser
add_subdirectory(external_libs)
#include_directories(${PROJECT_SOURCE_DIR}/external_libs/yaml-cpp/include)

# add include directory
#include_directories(${PROJECT_SOURCE_DIR})
include_directories(include cuda)

# list of source files
file(GLOB SOURCES "src/*.cpp")
if(CUDA_FOUND)
  file(GLOB HEADERS "include/*.h" "cuda/*.cuh")
else()
  file(GLOB HEADERS "include/*.h")
endif()
file(GLOB SOURCES_EXT_XML "external_libs/tinyxml2/*.cpp")

# add the executable
add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${SOURCES_EXT_XML})

# set include path
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE
          ${PROJECT_SOURCE_DIR}/include
          ${PROJECT_SOURCE_DIR}/cuda
          ${PROJECT_SOURCE_DIR}/external_libs/tinyxml2)

# compile cuda code
if (CUDA_FOUND)
  include_directories(${CUDA_INCLUDE_DIRS})
  add_library(${CUDA_LIBRARY_NAME} STATIC ${SOURCES_CUDA} ${HEADERS_CUDA} )
  target_include_directories(${CUDA_LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/external_libs/yaml-cpp/include)
  target_include_directories(${CUDA_LIBRARY_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_include_directories(${CUDA_LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/external_libs/tinyxml2)
endif()


# link mpi
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC MPI::MPI_CXX)

# link yaml-app:
target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC yaml-cpp)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/external_libs/yaml-cpp/include)

# link openmp
#target_link_libraries(Grid PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)
#target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)

# link HDF5
if(HDF5_FOUND)
  target_link_libraries(${CMAKE_PROJECT_NAME} PUBLIC ${HDF5_LIBRARIES})
  target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

# link cuda
if (CUDA_FOUND)

  #set_target_properties(${CUDA_LIBRARY_NAME} PROPERTIES CUDA_ARCHITECTURES "35;50;72")
  set_target_properties(${CUDA_LIBRARY_NAME} PROPERTIES CUDA_ARCHITECTURES "61")
  set_property(TARGET ${CUDA_LIBRARY_NAME} PROPERTY CUDA_ARCHITECTURES 61)

  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${CUDA_LIBRARY_NAME})
  target_link_libraries(${CUDA_LIBRARY_NAME} PUBLIC MPI::MPI_CXX)
  target_link_libraries(${CUDA_LIBRARY_NAME} PUBLIC yaml-cpp)
  target_link_libraries(${CUDA_LIBRARY_NAME} PUBLIC ${HDF5_LIBRARIES})

endif()