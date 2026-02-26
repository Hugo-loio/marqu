cmake_minimum_required(VERSION 3.10)
project(MarQu)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add library
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(marqu STATIC ${sources})
target_include_directories(marqu PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Install directives
install(TARGETS marqu DESTINATION "$ENV{HOME}/.local/lib")
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
install(FILES ${headers} DESTINATION "$ENV{HOME}/.local/include/${PROJECT_NAME}")

# Uninstall directive
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/uninstall.cmake)
