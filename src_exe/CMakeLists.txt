cmake_minimum_required(VERSION 3.0)
# Create Project
# THIS HAS TO COME BEFORE THE PROJECT LINE
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

# THIS HAS TO COME BEFORE THE PROJECT LINE
project (astarSolver)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -Wall )


#CACHE VARIABLES
set(PRINT_FILE OFF CACHE BOOL "")
set(USE_OPENMP ON CACHE BOOL "")
set(LOGGING ON CACHE BOOL "")
if(${LOGGING})
    add_compile_definitions(VERBOSE=1)
endif()
if(${PRINT_FILE})
    add_compile_definitions(ONFILE=1)
endif()

if(${USE_OPENMP})
    add_compile_definitions(USE_OPENMP=1)
    # OpenMP
    find_package(OpenMP REQUIRED COMPONENTS C)
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP_C_LIBRARIES: ${OpenMP_C_LIBRARIES}")
    message(STATUS "OpenMP_C_INCLUDE_DIRS: ${OpenMP_C_INCLUDE_DIRS}")

endif()




# ADD SOURCE SUBDIRECTORY
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../src_headers ${CMAKE_CURRENT_BINARY_DIR}/../aux-build)

# COPY HEADERS IN TEST INCLUDE FOLDER
set(headers_directory ${CMAKE_CURRENT_BINARY_DIR}/../include)
foreach(header ${aux_headers})
    get_filename_component(headerDirectory ${header} DIRECTORY)
    get_filename_component(headerFileName ${header} NAME_WE)
    get_filename_component(headerExtension ${header} EXT)
    message(${header})
    if(${headerExtension} MATCHES ".hpp.in")
        set(headerExtension ".hpp")

        configure_file(${header} ${headers_directory}/${headerFileName}${headerExtension})
    elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
        message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
    endif()

endforeach()


set(${PROJECT_NAME}_include ${${PROJECT_NAME}_include} ${headers_directory})


# INCLUDE_DIRECTORIES
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src_headers)



# ADD TEST EXECUTABLE
add_executable(${PROJECT_NAME} main.cpp ${aux_sources})


if(${USE_OPENMP})
    target_include_directories(${PROJECT_NAME} PRIVATE ${OpenMP_C_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}  ${OpenMP_C_LIBRARIES})
    target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_C_FLAGS})
endif()

