cmake_minimum_required(VERSION 3.10)

# Build options
#option(BUILD_OPENGL_EXAMPLE "Build OpenGL example." OFF)
set(CMAKE_BUILD_TYPE Release)
#add_compile_options(-g)

# Create the project
set(PROJECT_NAME sph_sim )
set(EIGEN_DIR /home/evgenii/libs/eigen-3.2.10)
set(SDL2_DIR /home/evgenii/libs/SDL2-2.0.10/)
project(${PROJECT_NAME})

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

#find all files that need to be compiled
file(GLOB ${PROJECT_NAME}_src ${PROJECT_ROOT}/src/*.cpp)
file(GLOB ${PROJECT_NAME}_h ${PROJECT_ROOT}/include/*.h)

#link_directories(${SDL2_DIR}/build)

#might want to compile with c++11 option or even 14
add_executable(${PROJECT_NAME}.exe ${${PROJECT_NAME}_src} ${${PROJECT_NAME}_h})

target_include_directories(${PROJECT_NAME}.exe PUBLIC ${PROJECT_ROOT}/include/ ${EIGEN_DIR} )
target_link_libraries(${PROJECT_NAME}.exe SDL2 SDL2_image)

# file(GLOB test_src ${PROJECT_ROOT}/src_other/*.cpp)
# file(GLOB test_h ${PROJECT_ROOT}/include_other/*.h)

# add_executable(test.exe ${test_src} ${test_h})
# target_include_directories(test.exe PUBLIC ${PROJECT_ROOT}/include_other/ ${EIGEN_DIR} )
# target_link_libraries(test.exe SDL2 SDL2_image)
