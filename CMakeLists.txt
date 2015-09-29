PROJECT(bemquad)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# C++11 required to build this library
set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

# put all the source files in one variable
file(GLOB bemquad_sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# put all the header files in one variables
file(GLOB bemquad_headers "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

# set the include directory for the header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

# create the bemquad library
ADD_LIBRARY(bemquad SHARED ${bemquad_sources} ${bemquad_headers})

# setup testing
INCLUDE(CTest)
ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/tests")
ENABLE_TESTING()
