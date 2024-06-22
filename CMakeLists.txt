cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(UAMMDstructured LANGUAGES CXX CUDA)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

# By default, install is done at ~/bin
set(INSTALL_DIR "$ENV{HOME}/bin" CACHE PATH "Installation directory")

set(UAMMD_REPO_URL "https://github.com/PabloIbannez/UAMMD.git")
set(UAMMD_TAG "v2.x")

###############################################################
# Download UAMMD
include(ExternalProject)

ExternalProject_Add(
    UAMMD
    GIT_REPOSITORY ${UAMMD_REPO_URL}
    GIT_TAG ${UAMMD_TAG}
    PREFIX ${CMAKE_BINARY_DIR}/external
    UPDATE_DISCONNECTED 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
    GIT_SHALLOW 1
)

###############################################################
# Create UAMMDstructured folder

set(UAMMD_STRUCTURED_PATH "${CMAKE_BINARY_DIR}/UAMMDstructured")
file(MAKE_DIRECTORY ${UAMMD_STRUCTURED_PATH})
file(MAKE_DIRECTORY ${UAMMD_STRUCTURED_PATH}/uammd)
file(MAKE_DIRECTORY ${UAMMD_STRUCTURED_PATH}/extensions)
###############################################################
# Copy structured/ folder to UAMMDstructured/

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/structured ${UAMMD_STRUCTURED_PATH}/extensions/structured SYMBOLIC)
###############################################################
# File generation

set(USCM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/USCM/USCM.py")

execute_process(COMMAND sh -c "python ${USCM_PATH} --generate"
                RESULT_VARIABLE result)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "USCM --generate failed.")
endif()

###############################################################
# Copy external/src/UAMMD/src/ to UAMMDstructured/src
# (done at build time)
# If UAMMDstructured/uammd/uammd.cuh is not found, copy the UAMMD source files to UAMMDstructured/uammd
add_custom_command(
    OUTPUT ${UAMMD_STRUCTURED_PATH}/uammd/uammd.cuh
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_BINARY_DIR}/external/src/UAMMD/src ${UAMMD_STRUCTURED_PATH}/uammd
    DEPENDS UAMMD
)

add_custom_target(
    CopyUAMMD ALL
    DEPENDS ${UAMMD_STRUCTURED_PATH}/uammd/uammd.cuh
)

add_dependencies(CopyUAMMD UAMMD)
################################################################
# Create UAMMDlauncher

# Set paths for generated files
set(LAUNCHER_DIR "${CMAKE_BINARY_DIR}/launcher")
file(MAKE_DIRECTORY ${LAUNCHER_DIR})

# Generate launcher code
file(GENERATE OUTPUT ${LAUNCHER_DIR}/UAMMDlauncher.cu CONTENT "
#include \"UAMMDstructured.cuh\"

using namespace uammd::structured;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        uammd::System::log<uammd::System::CRITICAL>(\"No input file provided!\");
        return EXIT_FAILURE;
    }

    std::string inputFilePath = argv[1];
    startSelfStartingSimulation(inputFilePath);

    return EXIT_SUCCESS;
}
")
################################################################
# Include necessary directories

# Include UAMMD
include_directories(${UAMMD_STRUCTURED_PATH}/uammd/)
# Include UAMMD third party
include_directories(${UAMMD_STRUCTURED_PATH}/uammd/third_party/)
# Include UAMMD structured
include_directories(${UAMMD_STRUCTURED_PATH}/extensions/structured/include/)
# Include CUDA
include_directories(${CUDA_INCLUDE_DIRS})

################################################################
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
################################################################
# Debug mode

# TODO: Add debug mode

################################################################
# Compilation flags
# Set C++ Standard
set(CMAKE_CXX_STANDARD 14)

################################################################
# CUDA Compiler options

# Add definitions: UAMMD_EXTENSIONS, UAMMD_EXTENSIONS_PREAMBLE=${UAMMD_STRUCTURED_PATH}/extensions/structured/preamble.h, MAXLOGLEVEL=5
add_definitions(-DUAMMD_EXTENSIONS -DUAMMD_EXTENSIONS_PREAMBLE=${UAMMD_STRUCTURED_PATH}/extensions/structured/preamble.h -DMAXLOGLEVEL=5)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --expt-extended-lambda -O3 -Xcompiler=\"-O3 -fPIC\" -w")

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

################################################################
# Packages

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

################################################################
# Source files
set(SOURCE_FILES ${DEST_DIR_MAIN}/launcher/UAMMDlauncher.cu)

# Recursively find all .cpp or .cu files in the structured directory and subdirectories
set(UAMMD_STRUCTURED_SRC_PATH "${UAMMD_STRUCTURED_PATH}/extensions/structured/src")
file(GLOB_RECURSE STRUCTURED_CPP_FILES "${UAMMD_STRUCTURED_SRC_PATH}/*.cpp")
file(GLOB_RECURSE STRUCTURED_CU_FILES "${UAMMD_STRUCTURED_SRC_PATH}/*.cu")

# Check if any .cu files were found
if(NOT STRUCTURED_CPP_FILES AND NOT STRUCTURED_CU_FILES)
    message(FATAL_ERROR "No .cpp or .cu files found in the structured directory.")
endif()

# Create an object library for the structured .cpp and .cu files
add_library(structured_objects STATIC ${STRUCTURED_CPP_FILES} ${STRUCTURED_CU_FILES})
add_dependencies(structured_objects CopyUAMMD)

# Define the executable
add_executable(UAMMDlauncher ${LAUNCHER_DIR}/UAMMDlauncher.cu)
add_dependencies(UAMMDlauncher CopyUAMMD)

# Link libraries
target_link_libraries(UAMMDlauncher PRIVATE
    structured_objects 
    $<TARGET_OBJECTS:structured_objects>
    ${LIBRARIES})

# Specify that the main build depends on the external project
add_dependencies(UAMMDlauncher UAMMD)

################################################################
## Installation rules
install(TARGETS UAMMDlauncher
        RUNTIME DESTINATION ${INSTALL_DIR})
