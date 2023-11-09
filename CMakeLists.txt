cmake_minimum_required(VERSION 3.26)
project(TSP_CP)

set(CMAKE_CXX_STANDARD 20)
################
##  OR_TOOLS  ##
################

# Fetch OR-tools library and create the alias ortools::ortools.
include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(BUILD_DEPS ON)
set(BUILD_SAMPLES OFF)
set(BUILD_EXAMPLES OFF)
set(ABSL_PROPAGATE_CXX_STD ON)
set(CPLEX_DIR /opt/ibm/ILOG/CPLEX_Studio2211)
include_directories(${CPLEX_DIR}/cplex/include)
include_directories(${CPLEX_DIR}/concert/include)
include_directories(${CPLEX_DIR}/cpoptimizer/include)

link_directories(${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic)
link_directories(${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic)
link_directories(${CPLEX_DIR}/cpoptimizer/lib/x86-64_linux/static_pic)

FetchContent_Declare(
        or-tools
        GIT_REPOSITORY https://github.com/google/or-tools.git
        GIT_TAG master

)
FetchContent_MakeAvailable(or-tools)
#set (HEADER_FILES
#        include/instance.h
#        include/cp_solver.h)

set (SOURCE_FILES
#        src/main.cpp src/instance.cpp
#        src/cp_solver.cpp)
        main.cpp
)
add_executable(TSP_CP ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(TSP_CP
        ortools::ortools
        ${CPLEX_DIR}/cpoptimizer/lib/x86-64_linux/static_pic/libcp.a
        ${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic/libconcert.a
        cp
        ilocplex
        ${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic/libcplex.a
        m
        pthread
        dl
)
