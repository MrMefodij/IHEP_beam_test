cmake_minimum_required(VERSION 3.21)
project(CubesLY)
set(CMAKE_CXX_STANDARD 17)

INCLUDE(${CMAKE_BINARY_DIR}/FindROOT.cmake)

set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++17")

INCLUDE_DIRECTORIES(. ${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(. ${CMAKE_CURRENT_SOURCE_DIR}/hitStructure/)
LINK_DIRECTORIES(     ${CMAKE_CURRENT_SOURCE_DIR}/lib)

if (APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    cmake_policy(SET CMP0042 NEW)
    cmake_policy(SET CMP0144 NEW)
    cmake_policy(SET CMP0153 NEW)
endif(APPLE)

set(CMAKE_PREFIX_PATH ${ROOT_INCLUDE_DIRS} ${ROOT_LIBRARY_DIR})

if    (ROOT_FOUND)
    INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
    LINK_DIRECTORIES   (${ROOT_LIBRARY_DIR})
endif (ROOT_FOUND)


message(STATUS "C flags      : " ${CMAKE_C_FLAGS})
message(STATUS "CXX flags    : " ${CMAKE_CXX_FLAGS})
message(STATUS "Linker flags : " ${CMAKE_EXE_LINKER_FLAGS}\n)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/hitStructure)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)

#install(FILES ${CMAKE_BINARY_DIR}/hitStructure/libPrototypeHitEventDict_rdict.pcm
#        ${CMAKE_BINARY_DIR}/hitStructure/libPrototypeHitEventDict.rootmap
#        ${CMAKE_BINARY_DIR}/hitStructure/PrototypeHitEventDict.cxx
#        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)