message(STATUS "examples ... ")

ADD_EXECUTABLE(eventStructure                                   SFGD_event_structure.cpp)
TARGET_LINK_LIBRARIES(eventStructure                            Gatchina_src PrototypeHitStructure)

ADD_EXECUTABLE(cubesLocation                                    SFGD_cubes_location.cpp)
TARGET_LINK_LIBRARIES(cubesLocation                             Gatchina_src PrototypeHitStructure)

ADD_EXECUTABLE(cubesResponse                                    SFGD_cubes_ly.cpp)
TARGET_LINK_LIBRARIES(cubesResponse                             Gatchina_src PrototypeHitStructure)

ADD_EXECUTABLE(mppcCalibration                                  SFGD_MPPC_calibration.cpp)
TARGET_LINK_LIBRARIES(mppcCalibration                           Gatchina_src PrototypeHitStructure)

ADD_EXECUTABLE(XTalk                                            SFGD_XTalk.cpp)
TARGET_LINK_LIBRARIES(XTalk                                     Gatchina_src PrototypeHitStructure)