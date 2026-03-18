cmake_minimum_required(VERSION 3.10)

project(MyGraphProject C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(external/GraphBLAS)
add_subdirectory(external/LAGraph)

add_executable(algoritms_bfs algoritms_bfs.c)
add_executable(main main.c)

target_link_libraries(algoritms_bfs PRIVATE GraphBLAS::GraphBLAS LAGraph::LAGraph)
target_link_libraries(main PRIVATE GraphBLAS::GraphBLAS LAGraph::LAGraph)

