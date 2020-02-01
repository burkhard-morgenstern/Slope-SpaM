cmake_minimum_required(VERSION 3.8.2)

project(slope-spam)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_subdirectory(third_party/fmt)

find_package(Threads)

include_directories(
    ${CMAKE_SOURCE_DIR}/third_party
    ${CMAKE_SOURCE_DIR}/third_party/fmt/include
)

add_library(slope-spam-common OBJECT
    src/spam/distance_matrix.cpp
    src/spam/pattern.cpp
    src/spam/sequence.cpp
    src/spam/wordlist.cpp
    src/config.cpp
    src/filesystem.cpp
    src/string.cpp
)

add_executable(slope-spam
    $<TARGET_OBJECTS:slope-spam-common>
    src/main.cpp
)

target_link_libraries(slope-spam
    ${CMAKE_THREAD_LIBS_INIT}
    fmt::fmt
    stdc++fs
)

# add_executable(slope-spam-test
#     $<TARGET_OBJECTS:slope-spam-common>
#     src/config.test.cpp
#     src/main.test.cpp
#     src/string.test.cpp
#     src/spam/pattern.test.cpp
#     src/spam/wordlist.test.cpp
# )

# target_include_directories(slope-spam-test PRIVATE
#     ${CMAKE_SOURCE_DIR}/third_party/catch2/single_include
# )

# target_link_libraries(slope-spam-test
#     ${CMAKE_THREAD_LIBS_INIT}
#     fmt::fmt
#     stdc++fs
# )

# add_test(slope-spam-test slope-spam-test)

# add_executable(generator
#     src/generator.cpp
# )

# add_executable(simulate
#     $<TARGET_OBJECTS:slope-spam-common>
#     src/simulate.cpp
#     src/generator.cpp
# )

# target_compile_definitions(simulate
#     PRIVATE -DGENERATOR_NO_MAIN
# )

# target_link_libraries(simulate
#     ${CMAKE_THREAD_LIBS_INIT}
#     fmt::fmt
#     stdc++fs
# )
