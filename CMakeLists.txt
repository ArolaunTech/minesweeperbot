cmake_minimum_required(VERSION 3.25)

project(configrender LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

set(
	sweeper_sources
	"src/main.cpp"
	"src/board.cpp"
	"src/board.h"
	"src/game.cpp"
	"src/game.h"
	"src/solver.cpp"
	"src/solver.h"
)

add_subdirectory("glfw-3.4")

add_library (glad STATIC "glad/src/glad.c")
target_include_directories(glad PUBLIC "glad/include")

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-Wall temp)

if(temp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

add_executable(sweepersolver ${sweeper_sources} ${PROJECT_PATH})
target_link_libraries(sweepersolver glfw glad)

target_include_directories(sweepersolver PUBLIC "glad/include")