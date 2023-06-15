add_executable(main main.cpp)
target_compile_options(main PRIVATE -Wall -Wextra -pedantic )
target_link_libraries(main  ${third_library})

