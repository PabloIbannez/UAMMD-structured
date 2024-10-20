cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(UAMMDstructured LANGUAGES CXX CUDA)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

# By default, install is done at ~/bin
set(INSTALL_DIR "$ENV{HOME}/bin" CACHE PATH "Installation directory")

# Add option for Python package installation
option(INSTALL_PYTHON_PACKAGE "Install Python package during make install" OFF)

# Add option for Python wrapper (only available if INSTALL_PYTHON_PACKAGE is ON)
option(BUILD_PYTHON_WRAPPER "Build Python wrapper for UAMMDlauncher" OFF)

# Add option for building .deb package
option(BUILD_DEB_PACKAGE "Build .deb package" OFF)

# Add option for setting ccbin for nvcc
set(CUDA_HOST_COMPILER "DEFAULT" CACHE STRING "Path to the host compiler for CUDA (DEFAULT to use default compiler)")

# Ensure BUILD_PYTHON_WRAPPER is OFF if INSTALL_PYTHON_PACKAGE is OFF
if(BUILD_PYTHON_WRAPPER AND NOT INSTALL_PYTHON_PACKAGE)
    message(FATAL_ERROR "BUILD_PYTHON_WRAPPER is ON but INSTALL_PYTHON_PACKAGE is OFF. "
                        "Please enable INSTALL_PYTHON_PACKAGE to build the Python wrapper.")
endif()

# BUILD_DEB_PACKAGE is only available if BUILD_PYTHON_WRAPPER and INSTALL_PYTHON_PACKAGE are OFF
if(BUILD_DEB_PACKAGE AND (BUILD_PYTHON_WRAPPER OR INSTALL_PYTHON_PACKAGE))
    message(FATAL_ERROR "BUILD_DEB_PACKAGE is ON but BUILD_PYTHON_WRAPPER or INSTALL_PYTHON_PACKAGE are ON. "
                        "Please disable BUILD_PYTHON_WRAPPER and INSTALL_PYTHON_PACKAGE to build the .deb package.")
endif()

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
set(UAMMD_COMPUTABLES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/structured/Computables.json")
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
# Set the precision of UAMMD, SINGLE or DOUBLE
set(UAMMD_PRECISION "SINGLE" CACHE STRING "Choose the precision of UAMMD, SINGLE or DOUBLE")

# If UAMMD_PRECISION is not SINGLE or DOUBLE, set it to SINGLE
if(NOT UAMMD_PRECISION STREQUAL "SINGLE" AND NOT UAMMD_PRECISION STREQUAL "DOUBLE")
  set(UAMMD_PRECISION "SINGLE")
  message(WARNING "UAMMD_PRECISION is not SINGLE or DOUBLE, set it to SINGLE")
endif()

################################################################
# Set Compiler Standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

################################################################
# CUDA Compiler options

# Set CMAKE_CUDA_HOST_COMPILER if CUDA_HOST_COMPILER is provided and not DEFAULT
if(NOT "${CUDA_HOST_COMPILER}" STREQUAL "DEFAULT")
    set(CMAKE_CUDA_HOST_COMPILER "${CUDA_HOST_COMPILER}")
endif()

################################################################
# Packages

# Check if MKLROOT is defined
if(DEFINED ENV{MKLROOT})
  message(STATUS "Detected MKLROOT: $ENV{MKLROOT}")
  set(BLA_VENDOR Intel10_64lp_seq)
  set(USE_MKL 1)
else()
  message(STATUS "MKLROOT not found. Using OpenBLAS.")
endif()

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

# Link libraries
set(LIBRARIES cufft cublas curand cusolver cusparse stdc++fs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(INSTALL_PYTHON_PACKAGE)
    find_package(Python COMPONENTS Interpreter REQUIRED)
endif()

if(BUILD_PYTHON_WRAPPER)
    # Find Python
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    # Find pybind11
    find_package(pybind11 REQUIRED)
endif()

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
    string(APPEND properties "((${cpp_type})(${name})(${c_type}))\\\n")
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

#define ${SECTION_UPPER}_PROPERTIES \\
${PROPERTIES}
")

  # Add a custom target for generating the header file
  add_custom_target(
    Generate${SECTION}Data ALL
    DEPENDS ${HEADER_FILE}
  )
endforeach()

# Generate computables

# Function to generate Computables.cuh
function(generate_computables_header json_file output_file)
  file(READ ${json_file} json_content)

  # Parse Computables
  string(JSON computables GET ${json_content} Computables)
  set(computables_list "")
  string(JSON computables_length LENGTH ${computables})
  math(EXPR computables_length "${computables_length} - 1")
  foreach(index RANGE ${computables_length})
    string(JSON item GET ${computables} ${index})
    string(JSON name GET ${item} 0)
    string(JSON type GET ${item} 1)
    string(APPEND computables_list "((${name})(${type}))\\\n")
  endforeach()

  # Remove the last newline and backslash
  string(REGEX REPLACE "\\\\n$" "" computables_list "${computables_list}")

  # Parse ComputablesCombinations
  string(JSON combinations GET ${json_content} ComputablesCombinations)
  set(combinations_list "")
  string(JSON combinations_length LENGTH ${combinations})
  math(EXPR combinations_length "${combinations_length} - 1")
  foreach(index RANGE ${combinations_length})
    string(JSON item GET ${combinations} ${index})
    string(JSON item_length LENGTH ${item})
    set(combination "")
    math(EXPR item_length "${item_length} - 1")

    foreach(item_index RANGE ${item_length})
      string(JSON element GET ${item} ${item_index})
      string(APPEND combination "(${element})")
    endforeach()

    string(APPEND combinations_list "(${combination})\\\n")
  endforeach()

  # Remove the last newline and backslash
  string(REGEX REPLACE "\\\\n$" "" combinations_list "${combinations_list}")

  #Write the header file
  file(WRITE ${output_file}
"#pragma once

#define COMPUTABLES \\
${computables_list}

#define COMPUTABLES_COMBINATIONS \\
${combinations_list}
")

endfunction()

set(COMPUTABLES_HEADER_FILE "${UAMMD_STRUCTURED_INCLUDE_PATH}/Definitions/Computables.cuh")

generate_computables_header(${UAMMD_COMPUTABLES_PATH} ${COMPUTABLES_HEADER_FILE})

# Add a custom target for generating the header file
add_custom_target(
  GenerateComputablesData ALL
  DEPENDS ${COMPUTABLES_HEADER_FILE}
)

################################################################
set(UAMMD_STRUCTURED_SRC_PATH "${UAMMD_STRUCTURED_PATH}/extensions/structured/src")
# Add backbone files to STRUCTURED_BACKBONE_FILES

list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Input/Input.cu")
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Input/InputEntryManager.cu")
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
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Interactor/ManyBodyNonBonded/ManyBodyNonBondedInteractor.cu")
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
list(APPEND STRUCTURED_BACKBONE_FILES "${UAMMD_STRUCTURED_SRC_PATH}/Utils/String/StringUtils.cu")
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
          # message(STATUS "Added component ${component_name} at src${current_path}/${component_type}/${component_file}")
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

# Create an object library for the structured .cpp and .cu files
add_library(structured_objects STATIC ${STRUCTURED_BACKBONE_FILES} ${STRUCTURED_COMPONENTS})
add_dependencies(structured_objects CopyUAMMD)

# Add necessary definitions to structured_objects
target_compile_definitions(structured_objects PRIVATE
    UAMMD_EXTENSIONS
    UAMMD_EXTENSIONS_PREAMBLE="${UAMMD_STRUCTURED_PATH}/extensions/structured/preamble.h"
    MAXLOGLEVEL=5
)

if(UAMMD_PRECISION STREQUAL "DOUBLE")
    target_compile_definitions(structured_objects PRIVATE DOUBLE_PRECISION)
endif()

if(USE_MKL)
    target_compile_definitions(structured_objects PRIVATE USE_MKL)
endif()

# Set include directories for structured_objects
target_include_directories(structured_objects PRIVATE
    ${UAMMD_STRUCTURED_PATH}/uammd/
    ${UAMMD_STRUCTURED_PATH}/uammd/third_party/
    ${UAMMD_STRUCTURED_PATH}/extensions/structured/include/
    ${UAMMD_STRUCTURED_PATH}/extensions/structured/src/
    ${CUDA_INCLUDE_DIRS}
)

# Set compile options for structured_objects
target_compile_options(structured_objects PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --expt-extended-lambda -O3 -Xcompiler=-O3,-fPIC -w>
)

# Set properties for structured_objects
set_target_properties(structured_objects PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
    CUDA_STANDARD 17
)

# Define the executable
add_executable(UAMMDlauncher ${CMAKE_CURRENT_SOURCE_DIR}/main/UAMMDlauncher.cu)
add_dependencies(UAMMDlauncher CopyUAMMD)

if(UAMMD_PRECISION STREQUAL "DOUBLE")
    target_compile_definitions(UAMMDlauncher PRIVATE DOUBLE_PRECISION)
endif()

if(USE_MKL)
    target_compile_definitions(UAMMDlauncher PRIVATE USE_MKL)
endif()

# Set include directories for UAMMDlauncher
target_include_directories(UAMMDlauncher PRIVATE
    ${UAMMD_STRUCTURED_PATH}/uammd/
    ${UAMMD_STRUCTURED_PATH}/uammd/third_party/
    ${UAMMD_STRUCTURED_PATH}/extensions/structured/include/
    ${UAMMD_STRUCTURED_PATH}/extensions/structured/src/
    ${CUDA_INCLUDE_DIRS}
)

# Set compile options for UAMMDlauncher
target_compile_options(UAMMDlauncher PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --expt-extended-lambda -O3 -Xcompiler=-O3,-fPIC -w>
)

# Set definitions for UAMMDlauncher
target_compile_definitions(UAMMDlauncher PRIVATE
    UAMMD_EXTENSIONS
    UAMMD_EXTENSIONS_PREAMBLE="${UAMMD_STRUCTURED_PATH}/extensions/structured/preamble.h"
    MAXLOGLEVEL=5
)

# Set properties for UAMMDlauncher
set_target_properties(UAMMDlauncher PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
    CUDA_STANDARD 17
)

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

if(BUILD_PYTHON_WRAPPER)
    # Add the Python module target
    pybind11_add_module(pyUAMMDlauncher ${CMAKE_CURRENT_SOURCE_DIR}/pyUAMMD/utils/launcher/UAMMDlauncher.cu)

    # Link libraries
    target_link_libraries(pyUAMMDlauncher PRIVATE
                          structured_objects
                          $<TARGET_OBJECTS:structured_objects>
                          ${LIBRARIES})

    # Set include directories for pyUAMMDlauncher
    target_include_directories(pyUAMMDlauncher PRIVATE
        ${UAMMD_STRUCTURED_PATH}/uammd/
        ${UAMMD_STRUCTURED_PATH}/uammd/third_party/
        ${UAMMD_STRUCTURED_PATH}/extensions/structured/include/
        ${UAMMD_STRUCTURED_PATH}/extensions/structured/src/
        ${CUDA_INCLUDE_DIRS}
    )

    # Set compile options for pyUAMMDlauncher
    target_compile_options(pyUAMMDlauncher PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr --expt-extended-lambda -O3 -Xcompiler=-O3,-fPIC -w>
    )

    # Set definitions for pyUAMMDlauncher
    target_compile_definitions(pyUAMMDlauncher PRIVATE
        UAMMD_EXTENSIONS
        UAMMD_EXTENSIONS_PREAMBLE="${UAMMD_STRUCTURED_PATH}/extensions/structured/preamble.h"
        MAXLOGLEVEL=5
    )

    if(UAMMD_PRECISION STREQUAL "DOUBLE")
        target_compile_definitions(pyUAMMDlauncher PRIVATE DOUBLE_PRECISION)
    endif()

    if(USE_MKL)
        target_compile_definitions(pyUAMMDlauncher PRIVATE USE_MKL)
    endif()

    # Set properties for pyUAMMDlauncher
    set_target_properties(pyUAMMDlauncher PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_STANDARD 17
        CUDA_STANDARD 17
    )

    # Add dependencies
    add_dependencies(pyUAMMDlauncher CopyUAMMD)
    foreach(SECTION ${SECTIONS})
        add_dependencies(pyUAMMDlauncher Generate${SECTION}Data)
    endforeach()
endif()

################################################################
# Build .deb package
if(BUILD_DEB_PACKAGE)
    include(InstallRequiredSystemLibraries)

    set(CPACK_PACKAGE_NAME "uammdstructured")
    set(CPACK_PACKAGE_VERSION "1.0.0")
    set(CPACK_PACKAGE_CONTACT "Pablo Ibanez Freire <p.ibanez.fre@gmail.com>")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "UAMMDstructured")
    set(CPACK_PACKAGE_VENDOR "Universidad Autonoma de Madrid")

    set(CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")

    include(CPack)

    install(TARGETS UAMMDlauncher
            RUNTIME DESTINATION bin)
endif()

if(NOT BUILD_DEB_PACKAGE)

    install(TARGETS UAMMDlauncher
            RUNTIME DESTINATION ${INSTALL_DIR})

    if(INSTALL_PYTHON_PACKAGE)
        add_custom_target(
            install_python_package
            COMMAND ${CMAKE_COMMAND} -E echo "Installing Python package..."
            COMMAND ${Python_EXECUTABLE} -m pip install ${CMAKE_SOURCE_DIR}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )

        # Make sure the Python package is installed after the main target and Python wrapper (if built)
        add_dependencies(install_python_package UAMMDlauncher)
        if(BUILD_PYTHON_WRAPPER)
            add_dependencies(install_python_package pyUAMMDlauncher)
        endif()

        # Add the Python package installation to the default install target
        install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target install_python_package)")
        if(BUILD_PYTHON_WRAPPER)
            # Installation rules for Python module
            install(TARGETS pyUAMMDlauncher
                    LIBRARY DESTINATION ${Python_SITELIB}/pyUAMMD/utils/launcher)
        endif()
    endif()
endif()

