enable_testing()

add_executable(YaoGraphTest YaoGraphTest.cpp)
target_link_libraries(YaoGraphTest GTest::gtest_main)

if (WITH_CAIRO)
    target_link_libraries(YaoGraphTest ${Cairomm_LIBRARIES})
endif ()

if (WITH_CGAL)
    target_link_libraries(YaoGraphTest ${CGAL_LIBRARIES})
endif ()

add_executable(PriorityQueueTest PriorityQueueTest.cpp)
target_link_libraries(PriorityQueueTest GTest::gtest_main)

add_executable(ListIndexTreeTest ListIndexTreeTest.cpp)
target_link_libraries(ListIndexTreeTest GTest::gtest_main)

add_executable(GeneratorTest GeneratorTest.cpp)
target_link_libraries(GeneratorTest GTest::gtest_main)