cmake_minimum_required(VERSION 3.10)

project(Math5600_Final)

add_executable(gaussian gaussian_elimination.cpp)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(gaussian PRIVATE Threads::Threads)

add_executable(mod_newton mod_newton.cpp)

add_executable(third_derivative third_derivative.cpp)