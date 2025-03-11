cmake_minimum_required(VERSION 3.9)

# Import local CMakeModules folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

project(hybridgwais C CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED 17)

#############################################################################
# C++ source files and target definition
#############################################################################
set(cpp_sources
    main.cpp
    version.cpp
    Args.cpp
    hybridsys/Buffer.h
    hybridsys/BufferAllocator.h
    hybridsys/BufferFactory.h
    hybridsys/Device.h
    hybridsys/DeviceCategory.cpp
    hybridsys/FPGA.cpp
    hybridsys/GPU.cpp
    hybridsys/Hybridsys.cpp
    hybridsys/Spinlock.h
    hybridsys/ThreadUtils.cpp
    ThreadPool.h
    MinMaxHeap.h
    StatusFile.h
    FPGAConfigurationGWAIS.cpp
    FPGAHandler.cpp
    GPUHandler.h
    GPUIDCreator.cpp
    ResultHandler.h
    SNPDB.cpp
    PlinkParser.cpp
    utils.h
    Method.cpp
    ResultView.h
    Progress.cpp
    Processor.cpp
    CPUProcessor.cpp
    # GPU
    GPUKernels.cu
    GPUEngine.cu
)


#############################################################################
# Set up Git version extraction
#############################################################################
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.template
    ${CMAKE_SOURCE_DIR}/version.h
    )

#include_directories("${CMAKE_BINARY_DIR}/generated")


#############################################################################
# Configure libraries
#############################################################################

# NVIDIA
option(GPU_SUPPORT "Use GPU Support if possible" ON)
cmake_policy(SET CMP0077 NEW)

if(GPU_SUPPORT)
    set(CMAKE_CUDA_ARCHITECTURES "60;70;80;86;90")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_70;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70;-gencode=arch=compute_80,code=sm_80;-gencode=arch=compute_86,code=sm_86")
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -O0 -g -G")
    set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O3 -lineinfo")

    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    add_executable(${PROJECT_NAME} ${cpp_sources} ${cuda_sources})
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_CUDA_GPU _FORCE_INLINES)
    target_link_libraries(${PROJECT_NAME} CUDA::nvml)

    set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_PROPAGATE_HOST_FLAGS OFF)
else()
    message("CUDA or NVML not found or deactivated. Compiling without GPU support!")
    set_source_files_properties(${cuda_sources} PROPERTIES LANGUAGE CXX)
    add_definitions(-x c++)
    add_executable(${PROJECT_NAME} ${cpp_sources} ${cuda_sources})
endif()

target_link_libraries(${PROJECT_NAME} rt)

# pthreads
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

# Zlib
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

# Alpha Data API
option(FPGA_SUPPORT "Use FPGA Support if possible" ON)
find_package(ADMXRC3) # optional!
if(ADMXRC3_FOUND AND FPGA_SUPPORT)
    target_link_libraries(${PROJECT_NAME} ${ADMXRC3_LIBRARY})
    include_directories(${ADMXRC3_INCLUDE_DIR})
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_AD_FPGA)
    message("Found ADMXRC3. Enabling FPGA support!")
else()
    message("Alpha Data API not found or deactivated. Compiling without FPGA support!")
endif()

# Boost
find_package(Boost 1.65.1 REQUIRED
    COMPONENTS
    filesystem
    program_options
)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

# Intel TBB
#find_package(TBB REQUIRED)
#include_directories(${TBB_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} ${TBB_LIBRARIES})
target_link_libraries(${PROJECT_NAME} tbb)
#target_compile_definitions(${PROJECT_NAME} PUBLIC ${TBB_DEFINITIONS} TBB_USE_EXCEPTIONS)

## Hybrid system API
#include(ExternalProject)
#ExternalProject_Add(HybridSys
#    GIT_REPOSITORY "git@github.com:lwienbrandt/hybridcomputer_api.git"
#    GIT_TAG "master"
#    GIT_SHALLOW true
#    SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/hybridsys-api"
#    BINARY_DIR "${CMAKE_BINARY_DIR}/deps/hybridsys-api"
#    UPDATE_COMMAND ""
#    INSTALL_COMMAND ""
#    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Release"
#)
#include_directories(${CMAKE_SOURCE_DIR}/deps/hybridsys-api)
#add_library(HybridSysLib SHARED IMPORTED)
#set_target_properties(HybridSysLib PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/deps/hybridsys-api/libHybridSys.a)

#add_dependencies(${PROJECT_NAME} HybridSys)
#target_link_libraries(${PROJECT_NAME} HybridSysLib)

# Common C++ libs project
#ExternalProject_Add(CommonCpp
#    GIT_REPOSITORY "git@github.com:lwienbrandt/cpp-common.git"
#    GIT_TAG "master"
#    GIT_SHALLOW true
#    SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cpp-common"
#    BINARY_DIR "${CMAKE_BINARY_DIR}/deps/cpp-common"
#    UPDATE_COMMAND ""
#    INSTALL_COMMAND ""
#    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Release"
#    )
#include_directories(${CMAKE_SOURCE_DIR}/deps/cpp-common)
#add_dependencies(${PROJECT_NAME} CommonCpp)
#add_library(CommonCppLib SHARED IMPORTED)
#set_target_properties(CommonCppLib PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/deps/cpp-common/libcpp-common.a)

#############################################################################
# CMake options and variables
#############################################################################

set(OPT_SCORE_TYPE	double	CACHE STRING "Score precision (double or float)")
target_compile_definitions(${PROJECT_NAME} PUBLIC SCORE_TYPE=${OPT_SCORE_TYPE})

# additional CXX options
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -ggdb -DNDEBUG")

set(CXX_WARNINGS "-Wall -Wextra -Wpedantic -Wshadow -Wno-deprecated-declarations -Wmissing-include-dirs -Wstrict-aliasing=2 -Wfloat-equal -Wcast-align")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_OPTIONS} ${cxx_compiler_opts} ${CXX_WARNINGS} -fopenmp -march=native")
