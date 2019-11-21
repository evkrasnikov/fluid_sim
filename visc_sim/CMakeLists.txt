cmake_minimum_required(VERSION 3.10)

# Build options
#option(BUILD_OPENGL_EXAMPLE "Build OpenGL example." OFF)
#set(CMAKE_BUILD_TYPE Debug)
add_compile_options(-g -O3)

# Create the project
set(PROJECT_NAME visc_sim )
set(EIGEN_DIR /home/evgenii/libs/eigen-3.2.10)
set(SDL2_DIR /home/evgenii/libs/SDL2-2.0.10/)
project(${PROJECT_NAME})
set(EXE_NAME ${PROJECT_NAME}.exe)

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

#find all files that need to be compiled
file(GLOB ${PROJECT_NAME}_src ${PROJECT_ROOT}/src/*.cpp)
file(GLOB ${PROJECT_NAME}_h ${PROJECT_ROOT}/include/*.h)

#link_directories(${SDL2_DIR}/build)

#might want to compile with c++11 option or even 14
# add_executable(${PROJECT_NAME}.exe ${${PROJECT_NAME}_src} ${${PROJECT_NAME}_h})

# target_include_directories(${PROJECT_NAME}.exe PUBLIC ${PROJECT_ROOT}/include/ ${EIGEN_DIR} )
# target_link_libraries(${PROJECT_NAME}.exe SDL2 SDL2_image)

file(GLOB ${PROJECT_NAME}_src ${PROJECT_ROOT}/src/*.cpp)
file(GLOB ${PROJECT_NAME}_h ${PROJECT_ROOT}/include/*.h)

add_executable(${EXE_NAME} ${${PROJECT_NAME}_src} )
target_include_directories(${EXE_NAME} PUBLIC ${PROJECT_ROOT}/include/ ${EIGEN_DIR} )
target_link_libraries(${EXE_NAME} SDL2 SDL2_image)
