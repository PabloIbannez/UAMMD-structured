cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(UAMMDstructured)
enable_language(CUDA)
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION OFF)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "all")
endif()

set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin CACHE PATH "Installation directory")

option(INSTALL_PYTHON_PACKAGE "Install Python package during make install" OFF)
option(BUILD_PYTHON_WRAPPER "Build Python wrapper for UAMMDlauncher" OFF)
option(BUILD_DEB_PACKAGE "Build .deb package" OFF)

if(BUILD_PYTHON_WRAPPER AND NOT INSTALL_PYTHON_PACKAGE)
    message(FATAL_ERROR "BUILD_PYTHON_WRAPPER is ON but INSTALL_PYTHON_PACKAGE is OFF. "
                        "Please enable INSTALL_PYTHON_PACKAGE to build the Python wrapper.")
endif()
if(BUILD_DEB_PACKAGE AND (BUILD_PYTHON_WRAPPER OR INSTALL_PYTHON_PACKAGE))
    message(FATAL_ERROR "BUILD_DEB_PACKAGE is ON but BUILD_PYTHON_WRAPPER or INSTALL_PYTHON_PACKAGE are ON. "
                        "Please disable BUILD_PYTHON_WRAPPER and INSTALL_PYTHON_PACKAGE to build the .deb package.")
endif()

set(UAMMD_REPO_URL "https://github.com/RaulPPelaez/UAMMD.git")
set(UAMMD_TAG "v2.8.0")
include(FetchContent)
FetchContent_Declare(
  uammd
  GIT_REPOSITORY ${UAMMD_REPO_URL}
  GIT_TAG        ${UAMMD_TAG}
  PATCH_COMMAND git apply ${PROJECT_SOURCE_DIR}/structured/patches/inline.patch
  UPDATE_DISCONNECTED TRUE
)
FetchContent_MakeAvailable(uammd)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

set(UAMMD_STRUCTURED_PATH "${PROJECT_SOURCE_DIR}")
set(UAMMD_STRUCTURED_INCLUDE_PATH "${UAMMD_STRUCTURED_PATH}/structured/include")
set(UAMMD_STRUCTURED_SRC_PATH "${UAMMD_STRUCTURED_PATH}/structured/src")
set(UAMMD_COMPONENTS_PATH "${PROJECT_SOURCE_DIR}/structured/Components.json")
set(UAMMD_DATA_PATH "${PROJECT_SOURCE_DIR}/structured/Data.json")
set(UAMMD_COMPUTABLES_PATH "${PROJECT_SOURCE_DIR}/structured/Computables.json")

set(UAMMD_PREAMBLE_PATH ${PROJECT_SOURCE_DIR}/structured/preamble.h)
set(COMPUTABLES_HEADER_FILE "${UAMMD_STRUCTURED_INCLUDE_PATH}/Definitions/Computables.cuh")

set(UAMMD_PRECISION "SINGLE" CACHE STRING "Choose the precision of UAMMD, SINGLE or DOUBLE")

if(NOT UAMMD_PRECISION STREQUAL "SINGLE" AND NOT UAMMD_PRECISION STREQUAL "DOUBLE")
  set(UAMMD_PRECISION "SINGLE")
  message(WARNING "UAMMD_PRECISION is not SINGLE or DOUBLE, set it to SINGLE")
endif()

include_directories(${UAMMD_STRUCTURED_INCLUDE_PATH})
include_directories(${UAMMD_STRUCTURED_SRC_PATH})
add_compile_definitions(
  UAMMD_EXTENSIONS
  UAMMD_EXTENSIONS_PREAMBLE=${UAMMD_PREAMBLE_PATH}
  MAXLOGLEVEL=5
)
# Ignore CUDA warnings about inline global kernels
set(CMAKE_CUDA_FLAGS "-diag-suppress=20050")

if(INSTALL_PYTHON_PACKAGE)
    find_package(Python COMPONENTS Interpreter REQUIRED)
endif()

if(BUILD_PYTHON_WRAPPER)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    find_package(pybind11 REQUIRED)
endif()

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


generate_computables_header(${UAMMD_COMPUTABLES_PATH} ${COMPUTABLES_HEADER_FILE})

add_custom_target(
  GenerateComputablesData ALL
  DEPENDS ${COMPUTABLES_HEADER_FILE}
)

file(GLOB_RECURSE STRUCTURED_BACKBONE_FILES
  "${UAMMD_STRUCTURED_SRC_PATH}/*.cu"
)
set(STRUCTURED_COMPONENTS "")

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
        else()
          message(FATAL_ERROR "Error parsing JSON: ${value} is not a valid file")
        endif()
      endif()
    endif()
  endforeach()
  set(STRUCTURED_COMPONENTS ${STRUCTURED_COMPONENTS} PARENT_SCOPE)
endfunction()

parse_json_and_add_files(${UAMMD_COMPONENTS_PATH} ${UAMMD_STRUCTURED_SRC_PATH})
list(REMOVE_DUPLICATES STRUCTURED_COMPONENTS)
add_library(structured_objects STATIC ${STRUCTURED_BACKBONE_FILES} ${STRUCTURED_COMPONENTS})
set_target_properties(structured_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)
uammd_setup_target(structured_objects)


if(UAMMD_PRECISION STREQUAL "DOUBLE")
    target_compile_definitions(structured_objects PRIVATE DOUBLE_PRECISION)
endif()

add_executable(UAMMDlauncher ${PROJECT_SOURCE_DIR}/main/UAMMDlauncher.cu)
uammd_setup_target(UAMMDlauncher)
if(UAMMD_PRECISION STREQUAL "DOUBLE")
    target_compile_definitions(UAMMDlauncher PRIVATE DOUBLE_PRECISION)
endif()

target_link_libraries(UAMMDlauncher PRIVATE
                      structured_objects
                      $<TARGET_OBJECTS:structured_objects>)


foreach(SECTION ${SECTIONS})
    add_dependencies(UAMMDlauncher Generate${SECTION}Data)
endforeach()

if(BUILD_PYTHON_WRAPPER)
    pybind11_add_module(pyUAMMDlauncher ${PROJECT_SOURCE_DIR}/pyUAMMD/utils/launcher/UAMMDlauncher.cu)
    uammd_setup_target(pyUAMMDlauncher)
    target_link_libraries(pyUAMMDlauncher PRIVATE
                          structured_objects
                          $<TARGET_OBJECTS:structured_objects>
			)

    if(UAMMD_PRECISION STREQUAL "DOUBLE")
        target_compile_definitions(pyUAMMDlauncher PRIVATE DOUBLE_PRECISION)
      endif()

    foreach(SECTION ${SECTIONS})
        add_dependencies(pyUAMMDlauncher Generate${SECTION}Data)
    endforeach()
endif()

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
