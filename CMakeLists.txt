cmake_minimum_required(VERSION 3.5.1)

project(slope-spam)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

find_package(Threads)

include_directories(
    ${CMAKE_SOURCE_DIR}/third_party
)

add_library(slope-spam-common OBJECT
    src/spam.cpp
)

add_executable(slope-spam
    $<TARGET_OBJECTS:slope-spam-common>
    src/config.cpp
    src/main.cpp
    src/math.cpp
)

target_link_libraries(slope-spam
    ${CMAKE_THREAD_LIBS_INIT}
    stdc++fs
)

add_executable(slope-spam-test
    $<TARGET_OBJECTS:slope-spam-common>
    src/main.test.cpp
    src/spam.test.cpp
)

target_include_directories(slope-spam-test PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party/catch2/single_include
)

target_link_libraries(slope-spam-test
    ${CMAKE_THREAD_LIBS_INIT}
    stdc++fs
)

add_test(slope-spam-test slope-spam-test)
