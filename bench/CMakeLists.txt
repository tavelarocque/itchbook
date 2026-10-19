add_executable(itchbook_bench bench_parse.cpp)
target_include_directories(itchbook_bench PRIVATE ${PROJECT_SOURCE_DIR}/tests)
target_link_libraries(itchbook_bench PRIVATE itchbook::core benchmark::benchmark)
itchbook_set_warnings(itchbook_bench)
