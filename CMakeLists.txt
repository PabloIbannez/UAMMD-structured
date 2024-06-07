cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(UAMMD_launcher LANGUAGES CXX CUDA)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

# Set the precision of UAMMD, SINGLE or DOUBLE
set(UAMMD_PRECISION "SINGLE" CACHE STRING "Choose the precision of UAMMD, SINGLE or DOUBLE")

# If UAMMD_PRECISION is not SINGLE or DOUBLE, set it to SINGLE
# If UAMMD_PRECISION is DOUBLE, set the flag -DDOUBLE_PRECISION
if(NOT UAMMD_PRECISION STREQUAL "SINGLE" AND NOT UAMMD_PRECISION STREQUAL "DOUBLE")
  set(UAMMD_PRECISION "SINGLE")
  message(WARNING "UAMMD_PRECISION is not SINGLE or DOUBLE, set it to SINGLE")
endif()
if(UAMMD_PRECISION STREQUAL "DOUBLE")
  add_definitions(-DDOUBLE_PRECISION)
endif()

# Check if UAMMD_PATH and UAMMD_STRUCTURED_PATH are defined
if(NOT DEFINED ENV{UAMMD_PATH})
  message(FATAL_ERROR "UAMMD_PATH is not defined")
endif()
if(NOT DEFINED ENV{UAMMD_STRUCTURED_PATH})
  message(FATAL_ERROR "UAMMD_STRUCTURED_PATH is not defined")
endif()

# Include directories
set(UAMMD_PATH $ENV{UAMMD_PATH} CACHE PATH "Path to UAMMD")
set(UAMMD_STRUCTURED_PATH $ENV{UAMMD_STRUCTURED_PATH} CACHE PATH "Path to UAMMD Structured")

include_directories(${UAMMD_PATH}/src/)
include_directories(${UAMMD_PATH}/src/third_party/)
include_directories(${UAMMD_STRUCTURED_PATH}/)
include_directories(${CUDA_INCLUDE_DIRS})

# Set C++ Standard
set(CMAKE_CXX_STANDARD 14)

# Find CUDA
find_package(CUDA REQUIRED)

# Check if MKLROOT is defined
if(DEFINED ENV{MKLROOT})
  set(BLA_VENDOR Intel10_64lp_seq)
  add_definitions(-DUSE_MKL)
endif()

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Link libraries
set(LIBRARIES cufft cublas curand cusolver cusparse stdc++fs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# CUDA Compiler options
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --expt-extended-lambda -O3 -DUAMMD_EXTENSIONS -DMAXLOGLEVEL=5 -Xcompiler=\"-O3 -fPIC\" -w")

# Specify the CUDA architecture(s)
set(CUDA_ARCHITECTURES 60 61 62 70 72 75 80 86 87 89 CACHE STRING "CUDA architectures to compile for")
# Filter out architectures that are not supported by the current CUDA version
execute_process(
    COMMAND nvcc --help
    OUTPUT_VARIABLE NVCC_OUTPUT
)

# Parse the NVCC_OUTPUT to extract supported architectures
# This will need to be adjusted based on the actual format of the nvcc output
string(REGEX MATCHALL "compute_[0-9]+" SUPPORTED_ARCHS "${NVCC_OUTPUT}")

# Convert the architecture names to the format used in CUDA_ARCHITECTURES
set(PROCESSED_ARCHS "")
foreach(ARCH IN LISTS SUPPORTED_ARCHS)
    string(REGEX REPLACE "compute_" "" ARCH_NUM "${ARCH}")
    list(APPEND PROCESSED_ARCHS "${ARCH_NUM}")
endforeach()

# Filter CUDA_ARCHITECTURES based on the supported architectures
set(FILTERED_CUDA_ARCHITECTURES "")
foreach(ARCH IN LISTS CUDA_ARCHITECTURES)
    if(ARCH IN_LIST PROCESSED_ARCHS)
        list(APPEND FILTERED_CUDA_ARCHITECTURES ${ARCH})
    else()
        message(STATUS "Architecture ${ARCH} is not supported by the current NVCC version and has been removed.")
    endif()
endforeach()

# Set the filtered architectures
set(CMAKE_CUDA_ARCHITECTURES ${FILTERED_CUDA_ARCHITECTURES})

# Libraries
set(LIBRARIES cufft cublas curand cusolver cusparse stdc++fs)

# Source files
set(SOURCE_FILES ${UAMMD_PATH}/launcher/UAMMDlauncher.cu)

# Recursively find all .cu files in the structured directory and subdirectories
file(GLOB_RECURSE STRUCTURED_CU_FILES ${UAMMD_PATH}/extensions/structured/src/*.cu)

# Check if any .cu files were found
if(NOT STRUCTURED_CU_FILES)
  message(FATAL_ERROR "No .cu files found in ${UAMMD_PATH}/extensions/structured/src/")
endif()

# Print all paths
message(STATUS "UAMMD_PATH: " ${UAMMD_PATH})
message(STATUS "UAMMD_STRUCTURED_PATH: " ${UAMMD_STRUCTURED_PATH})
message(STATUS "CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})
message(STATUS "UAMMD_STRUCTURED OBJECTS: " ${STRUCTURED_CU_FILES})

# Define the executable
add_executable(UAMMDlauncher ${SOURCE_FILES})

# Create an object library for the structured .cu files
cuda_add_library(structured_objects SHARED ${STRUCTURED_CU_FILES})

# Link libraries
target_link_libraries(UAMMDlauncher structured_objects ${LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Installation rules
install(TARGETS UAMMDlauncher
        RUNTIME DESTINATION ${UAMMD_PATH}/bin)

if(UNIX)
    # For Unix-like systems
    execute_process(COMMAND sh -c "USCM --generate"
                    RESULT_VARIABLE result)
elseif(WIN32)
    # For Windows
    execute_process(COMMAND cmd /c "USCM --generate"
                    RESULT_VARIABLE result)
endif()

if(NOT result EQUAL 0)
  message(FATAL_ERROR "USCM --generate failed.")
endif()

