cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(CPLEX)
include_directories(${CPLEX_INCLUDE_DIRS})


#set(CMAKE_CXX_FLAGS "-zstack-size=8684354560 -frounding-math -DIL_STD -g3 -pedantic -Wall -Wextra -Wshadow -fsanitize=address -fno-sanitize-recover=all -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS "-frounding-math -DIL_STD -O3")

find_package(CGAL  REQUIRED COMPONENTS Core)

find_package(Boost COMPONENTS thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

link_libraries(CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES} ${CPLEX_LIBRARIES})

add_executable(BC src/main.cpp src/Arrangement.cpp src/SEC.cpp src/Solver.cpp src/DisjointSetForest.cpp src/Timer.cpp src/Benders.cpp)

add_executable(BC-PP src/main.cpp src/Arrangement.cpp src/SEC.cpp src/Solver.cpp src/DisjointSetForest.cpp src/Timer.cpp src/Benders.cpp)
target_compile_definitions(BC-PP PUBLIC PREPROCESS)

add_executable(BC-TRI src/main.cpp src/Arrangement.cpp src/SEC.cpp src/Solver.cpp src/DisjointSetForest.cpp src/Timer.cpp src/Benders.cpp)
target_compile_definitions(BC-TRI PUBLIC SEPARATE_TRIANGLE_CUTS)

add_executable(BC-TRI-PP src/main.cpp src/Arrangement.cpp src/SEC.cpp src/Solver.cpp src/DisjointSetForest.cpp src/Timer.cpp src/Benders.cpp)
target_compile_definitions(BC-TRI-PP PUBLIC PREPROCESS)
target_compile_definitions(BC-TRI-PP PUBLIC SEPARATE_TRIANGLE_CUTS)

add_executable(BC-BD src/main.cpp src/Arrangement.cpp src/SEC.cpp src/Solver.cpp src/DisjointSetForest.cpp src/Timer.cpp src/Benders.cpp)
target_compile_definitions(BC-BD PUBLIC BENDERS)


