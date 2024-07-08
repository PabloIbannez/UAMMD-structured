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

set(UAMMD_COMPONENTS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/structured/Components.json")
set(UAMMD_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/structured/Data.json")
###############################################################
# Copy structured/ folder to UAMMDstructured/

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/structured ${UAMMD_STRUCTURED_PATH}/extensions/structured SYMBOLIC)

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
# Include necessary directories

# Include UAMMD
include_directories(${UAMMD_STRUCTURED_PATH}/uammd/)
# Include UAMMD third party
include_directories(${UAMMD_STRUCTURED_PATH}/uammd/third_party/)
# Include UAMMD structured
include_directories(${UAMMD_STRUCTURED_PATH}/extensions/structured/include/)
include_directories(${UAMMD_STRUCTURED_PATH}/extensions/structured/src/)
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
# Files generation

set(UAMMD_STRUCTURED_INCLUDE_PATH "${UAMMD_STRUCTURED_PATH}/extensions/structured/include")

# Function to parse JSON and generate property definitions
function(generate_property_definitions json_file section output_var)
  file(READ ${json_file} json_content)
  string(JSON section_data GET ${json_content} ${section})

  set(properties "")
  string(JSON array_length LENGTH ${section_data})
  math(EXPR array_length "${array_length} - 1")

  foreach(index RANGE ${array_length})
    string(JSON item GET ${section_data} ${index})
    string(JSON name GET ${item} 0)
    string(JSON cpp_type GET ${item} 1)
    string(JSON c_type GET ${item} 2)
    string(APPEND properties "((${cpp_type},${name},${c_type}))\\\n")
  endforeach()

  # Remove the last newline and backslash
  string(REGEX REPLACE "\\\\n$" "" properties "${properties}")

  set(${output_var} "${properties}" PARENT_SCOPE)
endfunction()

# List of sections to process
set(SECTIONS ParticleData State Units Fundamental Ensemble)

# Generate header files for each section
foreach(SECTION ${SECTIONS})
  string(TOUPPER ${SECTION} SECTION_UPPER)
  set(HEADER_FILE "${UAMMD_STRUCTURED_INCLUDE_PATH}/Definitions/${SECTION}.cuh")

  generate_property_definitions(${UAMMD_DATA_PATH} ${SECTION} PROPERTIES)

  file(WRITE ${HEADER_FILE}
"#pragma once

#define ${SECTION_UPPER}_PROPERTIES ${PROPERTIES}
")

  # Add a custom target for generating the header file
  add_custom_target(
    Generate${SECTION}Data ALL
    DEPENDS ${HEADER_FILE}
  )
endforeach()

################################################################
set(UAMMD_STRUCTURED_SRC_PATH "${UAMMD_STRUCTURED_PATH}/extensions/structured/src")
# Add backbone files to STRUCTURED_BACKBONE_FILES
# Add System/ExtendedSystem.cu
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/System/ExtendedSystem.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/GlobalData/GlobalData.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/GlobalData/Units/UnitsLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/GlobalData/Ensemble/EnsembleLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/GlobalData/Types/TypesLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/GlobalData/Fundamental/FundamentalLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/ParticleData/ExtendedParticleData.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/ParticleData/StateLoader.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/ParticleGroup/ParticleGroupUtils.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Interactor/InteractorLoader.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Interactor/Patches/PatchInteractorLoader.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Interactor/PatchyParticles/PatchyParticlesInteractor.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Integrator/IntegratorBase.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Integrator/IntegratorManager.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Integrator/IntegratorLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Integrator/IntegratorUtils.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/ForceFields/ForceFields.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Topology/Topology.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/SimulationStep/SimulationStep.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/SimulationStep/SimulationStepManager.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/SimulationStep/GenericSimulationStepLoader.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Simulation/Simulation.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/DataStructures/VerletConditionalListSet/VerletConditionalListSetUtils.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/DataStructures/VerletConditionalListSet/VerletConditionalListSetLoaders.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/DataStructures/GroupsList/GroupsList.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/DataStructures/ExclusionsList/ExclusionsList.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/Output/Output.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/Backup/Backup.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/Backup/BackupStep.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/Measures/MeasuresBasic.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/ParameterHandler/CheckDataConsistency.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/Containers/SetUtils.cu")

# Define STRUCTURED_COMPONENTS as a global variable
set(STRUCTURED_COMPONENTS "")

# Add this function to parse the JSON and populate STRUCTURED_COMPONENTS
function(parse_json_and_add_files json_file base_path)
  file(READ ${json_file} json_content)
  string(JSON json_object ERROR_VARIABLE json_error GET ${json_content})

  if(json_error)
    message(FATAL_ERROR "Error parsing JSON: ${json_error}")
  endif()

  parse_json_object("" ${json_object})

  set(STRUCTURED_COMPONENTS ${STRUCTURED_COMPONENTS} PARENT_SCOPE)
endfunction()

function(parse_json_object current_path object)
  string(JSON num_elements LENGTH ${object})
  math(EXPR num_elements "${num_elements} - 1")

  foreach(i RANGE ${num_elements})
    string(JSON key MEMBER ${object} ${i})
    string(JSON value GET ${object} ${key})

    if("${value}" MATCHES "^\\[")
      parse_json_array("${current_path}/${key}" ${value})
    elseif("${value}" MATCHES "^\\{")
      parse_json_object("${current_path}/${key}" ${value})
    endif()
  endforeach()

  set(STRUCTURED_COMPONENTS ${STRUCTURED_COMPONENTS} PARENT_SCOPE)
endfunction()

function(parse_json_array current_path array)
  string(JSON num_elements LENGTH ${array})
  math(EXPR num_elements "${num_elements} - 1")

  foreach(i RANGE ${num_elements})
    string(JSON value GET ${array} ${i})
    if("${value}" MATCHES "^\\[")
      string(JSON inner_array_length LENGTH ${value})
      if(inner_array_length GREATER 2)
        string(JSON component_type GET ${value} 0)
        string(JSON component_name GET ${value} 1)
        string(JSON component_file GET ${value} 2)
        if(component_file MATCHES "\\.cu$")
          set(file_path "${base_path}${current_path}/${component_type}/${component_file}")
          list(APPEND STRUCTURED_COMPONENTS ${file_path})
          set(STRUCTURED_COMPONENTS ${STRUCTURED_COMPONENTS} PARENT_SCOPE)
          message(STATUS "Added component ${component_name} at src${current_path}/${component_type}/${component_file}")
        else()
          message(FATAL_ERROR "Error parsing JSON: ${value} is not a valid file")
        endif()
      endif()
    endif()
  endforeach()
  set(STRUCTURED_COMPONENTS ${STRUCTURED_COMPONENTS} PARENT_SCOPE)
endfunction()

# Use the function to populate STRUCTURED_COMPONENTS
parse_json_and_add_files(${UAMMD_COMPONENTS_PATH} ${UAMMD_STRUCTURED_PATH}/extensions/structured/src/)

# Remove duplicates in STRUCTURED_COMPONENTS
list(REMOVE_DUPLICATES STRUCTURED_COMPONENTS)

# Print the components
#message(STATUS "Components files:")
#foreach(component IN LISTS STRUCTURED_COMPONENTS)
#  message(STATUS ${component})
#endforeach()

# Create an object library for the structured .cpp and .cu files
add_library(structured_objects STATIC ${STRUCTURED_BACKBONE_FILES} ${STRUCTURED_COMPONENTS})
add_dependencies(structured_objects CopyUAMMD)
foreach(SECTION ${SECTIONS})
  add_dependencies(structured_objects Generate${SECTION}Data)
endforeach()

# Define the executable
add_executable(UAMMDlauncher ${CMAKE_CURRENT_SOURCE_DIR}/main/UAMMDlauncher.cu)
add_dependencies(UAMMDlauncher CopyUAMMD)

# Link libraries
target_link_libraries(UAMMDlauncher PRIVATE
                      structured_objects
                      $<TARGET_OBJECTS:structured_objects>
                      ${LIBRARIES})

# Specify that the main build depends on the external project and generated files
add_dependencies(UAMMDlauncher CopyUAMMD)
foreach(SECTION ${SECTIONS})
  add_dependencies(UAMMDlauncher Generate${SECTION}Data)
endforeach()

################################################################
## Installation rules
install(TARGETS UAMMDlauncher
        RUNTIME DESTINATION ${INSTALL_DIR})
