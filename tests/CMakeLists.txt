ADD_EXECUTABLE(detest DblExpQuadTest.cpp)
ADD_EXECUTABLE(adaptiveTrapTest adaptiveTrap.cpp)
ADD_EXECUTABLE(adaptive2DCCTest adaptive2DCCTest.cpp)

TARGET_LINK_LIBRARIES(detest bemquad)
TARGET_LINK_LIBRARIES(adaptiveTrapTest bemquad)
TARGET_LINK_LIBRARIES(adaptive2DCCTest bemquad)

ADD_TEST(integrationTest1 detest)
ADD_TEST(integrationTest2 adaptiveTrapTest )
ADD_TEST(integrationTest3 adaptive2DCCTest )
