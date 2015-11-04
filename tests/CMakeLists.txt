ADD_EXECUTABLE(detest DblExpQuadTest.cpp)
ADD_EXECUTABLE(adaptiveTrapTest adaptiveTrap.cpp)
ADD_EXECUTABLE(adaptive2DCCTest adaptive2DCCTest.cpp)
ADD_EXECUTABLE(complexQuadTest complexQuadratureTest.cpp)

TARGET_LINK_LIBRARIES(detest bemquad)
TARGET_LINK_LIBRARIES(adaptiveTrapTest bemquad)
TARGET_LINK_LIBRARIES(adaptive2DCCTest bemquad)
TARGET_LINK_LIBRARIES(complexQuadTest bemquad)

ADD_TEST(integrationTest1 detest)
ADD_TEST(integrationTest2 adaptiveTrapTest )
ADD_TEST(integrationTest3 adaptive2DCCTest )
ADD_TEST(integrationTest4 complexQuadTest)
