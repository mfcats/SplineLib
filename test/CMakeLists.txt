set(TEST_SOURCES
    example_test.cc)

find_package(GTest REQUIRED)

add_executable(SplineLibTests ${TEST_SOURCES})
target_link_libraries(SplineLibTests GTest::Main)