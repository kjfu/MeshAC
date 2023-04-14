#pragma once
#include <string>
#include "mesh.h"
#include "SurfaceMesh.h"
namespace MeshAC{

    void generateMeshFromPoints(
        const std::string &input, 
        const std::string &output, 
        double size
    );

    void generateMeshFromEdgeDislocationPoints(
        const std::string &input, 
        const std::string &output, 
        double size
    );

    void generateMeshForAtomisticRegion(
        const std::vector<Vector3D> &points,
        Mesh &resultingMesh,
        SurfaceMesh &interfaceMesh
    );

    void generateMeshForContinuumRegion(
        SurfaceMesh &boundingBoxSurfaceMesh,
        SurfaceMesh &interfaceSurfaceMesh,
        Mesh &resultingMesh
    );

    void generateBoundingBoxSurfaceMesh(
        const Vector3D &lowerBound,
        const Vector3D &upperBound,
        double size,
        SurfaceMesh &resultingSurfaceMesh
    );

    void removeTopBottomTrianglesInSurfaceMesh(
        SurfaceMesh &surfaceMesh,
        double top,
        double bottom,
        std::vector<std::array<double, 2>> &topEdgeNodes,
        std::vector<std::array<double, 2>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges,
        std::vector<std::array<int, 2>> &bottomEdges
    );

    void generateBoundingBoxSurfaceMeshWithXYHoles(
        const Vector3D &lowerBound,
        const Vector3D &upperBound,
        double size,
        std::vector<std::array<double, 2>> &topEdgeNodes,
        std::vector<std::array<double, 2>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges, 
        std::vector<std::array<int, 2>> &bottomEdges,
        SurfaceMesh &resultingSurfaceMesh
    );




}