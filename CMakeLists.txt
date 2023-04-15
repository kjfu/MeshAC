cmake_minimum_required(VERSION 3.5)
project(MeshAC)

set(CMAKE_CXX_STANDARD 14)

include_directories(
    ./src
    ./external/tetgen
    ./external/kdtree
    ./external/AABB
    ./external/triangle
)

set(MESHAC_SRCS
    ./src/main.cpp
    ./src/mesh.cpp
    ./src/mesher3d_core.cpp
    ./src/mesher3d_io.cpp
    ./src/MeshRefiner.cpp
    ./src/MeshSwapper.cpp
    ./src/SpatialSearcher.cpp
    ./src/sphere.cpp
    ./src/SurfaceMesh.cpp
    ./src/PointSet.cpp
    ./src/CommandParser.cpp
    ./src/TetgenTool.cpp
    ./src/TriangleTool.cpp
    ./src/Core.cpp
    ./src/CoreIO.cpp
)

#tetgen
set(EXTERNAL_TETGEN_SRCS
    ./external/tetgen/tetgen.cxx
    ./external/tetgen/predicates.cxx
)

#triangle
set(EXTERNAL_TRIANGLE_SRCS
    ./external/triangle/triangle.c
)

#kdtree
set(EXTERNAL_KDTREE_SRCS
    ./external/kdtree/kdtree.c
)

# AABB
set(EXTERNAL_AABB_SRCS
    ./external/AABB/AABB.cc
)



add_executable(MeshAC 
    ${MESHAC_SRCS} 
    ${EXTERNAL_KDTREE_SRCS}
    ${EXTERNAL_TETGEN_SRCS}
    ${EXTERNAL_AABB_SRCS}
    ${EXTERNAL_TRIANGLE_SRCS}
)

set_target_properties(MeshAC PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)