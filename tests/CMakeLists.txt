add_executable(itchbook_tests
  test_messages.cpp
  test_stream.cpp
  test_properties.cpp
  test_fields.cpp
  test_book.cpp
  test_bench_support.cpp
  test_queue.cpp
)
target_include_directories(itchbook_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(itchbook_tests PRIVATE itchbook::core GTest::gtest_main rapidcheck rapidcheck_gtest Threads::Threads)
itchbook_set_warnings(itchbook_tests)

include(GoogleTest)
gtest_discover_tests(itchbook_tests DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 300)

# Cross-check: the C++ fingerprint and the independent Python reference must
# agree byte-for-byte on a synthetic stream that contains every message type.
find_package(Python3 3.10 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_test(NAME fingerprint_matches_python_reference
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/itchref/crosscheck.py
            --itch-stats $<TARGET_FILE:itch-stats> --synthetic 20000
  )
  add_test(NAME queue_study_end_to_end
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/itchref/queue_smoke.py
            --itch-queue $<TARGET_FILE:itch-queue>
  )
  add_test(NAME book_matches_python_reference
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/itchref/crosscheck_book.py
            --itch-book $<TARGET_FILE:itch-book> --synthetic 50000
  )
endif()
