cmake_minimum_required(VERSION 3.14)
project(algo)

set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

FIND_PATH(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
        PATHS ${CMAKE_SOURCE_DIR}/lib/lemon-1.3.1 ${CMAKE_SOURCE_DIR}/deps/lemon
        NO_DEFAULT_PATH
        DOC "Location of LEMON source as a CMAKE subproject")

IF(EXISTS ${LEMON_SOURCE_ROOT_DIR})
    ADD_SUBDIRECTORY(${LEMON_SOURCE_ROOT_DIR} deps/lemon)
    SET(LEMON_INCLUDE_DIRS
            ${LEMON_SOURCE_ROOT_DIR}
            ${CMAKE_BINARY_DIR}/deps/lemon
            )
    SET(LEMON_LIBRARIES lemon)
    UNSET(LEMON_ROOT_DIR CACHE)
    UNSET(LEMON_DIR CACHE)
    UNSET(LEMON_INCLUDE_DIR CACHE)
    UNSET(LEMON_LIBRARY CACHE)
ELSE()
    FIND_PACKAGE(LEMON QUIET NO_MODULE)
    FIND_PACKAGE(LEMON REQUIRED)
ENDIF()

## This line finds doxygen (for document creation)

FIND_PACKAGE(Doxygen)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_BINARY_DIR}
        ${LEMON_INCLUDE_DIRS}
)

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


add_subdirectory(lib/pybind11-master)

#pybind11_add_module(pywrap UniqueID.cpp UniqueID.h Specie.cpp Specie.h ContactNetwork.cpp ContactNetwork.h SSA.cpp SSA.h py_wrap.cpp)
add_executable(algo main.cpp UniqueID.cpp UniqueID.h Specie.cpp Specie.h ContactNetwork.cpp ContactNetwork.h SSA.cpp SSA.h NSA.cpp NSA.h)