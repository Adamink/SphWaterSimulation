cmake_minimum_required(VERSION 3.9)
project(SphWaterSimulation)

# set(CMAKE_CXX_FLAGS "-Wall")

# source files
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SphWaterSimulation")
file(GLOB SRCFILES ${SRC_DIR}/*.cc ${SRC_DIR}/*.cpp ${SRC_DIR}/*.c ${SRC_DIR}/*.hpp ${SRC_DIR}/*.h)
add_executable(${PROJECT_NAME} ${SRCFILES})

set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")

# freeglut
set(GLUT_DIR "${LIB_DIR}/freeglut")
set(FREEGLUT_BUILD_DEMOS OFF CACHE INTERNAL "Build FreeGLUT demos")
set(FREEGLUT_BUILD_STATIC_LIBS OFF CACHE INTERNAL " Build FreeGLUT static library")
add_subdirectory("${GLUT_DIR}")
target_link_libraries(${PROJECT_NAME} "freeglut" "GLU")
target_include_directories(${PROJECT_NAME} PRIVATE "${GLUT_DIR}/include")

# omp
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()