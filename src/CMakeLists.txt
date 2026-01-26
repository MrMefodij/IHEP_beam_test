message(STATUS "src ... ")

AUX_SOURCE_DIRECTORY(.                 SOURCE_FILES)

if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})

endif()

ADD_LIBRARY(Gatchina_src      SHARED            ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(Gatchina_src              ${Boost_LIBRARIES} ${ROOT_LIBRARIES} -lSpectrum  PrototypeHitStructure)