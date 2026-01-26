message(STATUS "hitStructure")
AUX_SOURCE_DIRECTORY(.                 SOURCE_FILES)
#message(STATUS "${SOURCE_FILES}")
INCLUDE_DIRECTORIES(. )

message(STATUS "Looking for ROOT")
find_package(ROOT REQUIRED COMPONENTS RIO)
message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}")
include_directories(${ROOT_INCLUDE_DIRS} .)

ROOT_GENERATE_DICTIONARY(PrototypeHitEventDict PrototypeHit.h LINKDEF LinkDef.h )

ADD_LIBRARY(PrototypeHitStructure      SHARED            ${SOURCE_FILES}  PrototypeHitEventDict.cxx LinkDef.h)
TARGET_LINK_LIBRARIES(PrototypeHitStructure     PUBLIC    ${ROOT_LIBRARIES})

#add_custom_command(
#        TARGET PrototypeHitStructure
#        POST_BUILD
#        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
#)