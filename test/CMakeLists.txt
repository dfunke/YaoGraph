find_package(GTest)

add_executable(GTest_run SweepLineTest.cpp ListIndexTreeTest.cpp)
target_link_libraries(GTest_run GTest::gtest_main)

if (WITH_CAIRO)
    target_link_libraries(GTest_run ${Cairomm_LIBRARIES})
endif ()

if (WITH_CGAL)
    target_link_libraries(GTest_run ${CGAL_LIBRARIES})
endif ()