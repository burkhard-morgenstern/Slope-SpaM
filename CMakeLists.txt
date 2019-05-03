cmake_minimum_required(VERSION 3.10.2)

project(slope-spam)

set(CMAKE_CXX_STANDARD 17)

find_package(Threads)

add_executable(slope-spam
    src/config.cpp
    src/main.cpp
    src/math.cpp
    src/spam.cpp
)

target_include_directories(slope-spam PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party
)

target_link_libraries(slope-spam
    ${CMAKE_THREAD_LIBS_INIT}
)
