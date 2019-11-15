cmake_minimum_required(VERSION 3.10)

# Build options
#option(BUILD_OPENGL_EXAMPLE "Build OpenGL example." OFF)

# Create the project
set(PROJECT_NAME sph_sim )
set(EIGEN_DIR /home/evgenii/libs/eigen-3.2.10)
project(${PROJECT_NAME})

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

#find all files that need to be compiled
file(GLOB ${PROJECT_NAME}_src ${PROJECT_ROOT}/src/*.cpp)
file(GLOB ${PROJECT_NAME}_h ${PROJECT_ROOT}/include/*.h)

#might want to compile with c++11 option or even 14
add_executable(${PROJECT_NAME}.exe ${${PROJECT_NAME}_src} ${${PROJECT_NAME}_h})

target_include_directories(${PROJECT_NAME}.exe PUBLIC ${PROJECT_ROOT}/include/ ${EIGEN_DIR})

