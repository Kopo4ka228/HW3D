cmake_minimum_required(VERSION 3.5)

project(triangles)

set(SOURCES
    funcs.cpp
    main.cc
)

set(HEADERS
    geom.hpp
)

add_executable(triangles ${SOURCES} ${HEADERS})
set_target_properties(
    triangles PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)