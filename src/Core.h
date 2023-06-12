/*
 * @Author: Kejie Fu
 * @Date: 2023-04-13 23:16:47
 * @LastEditTime: 2023-06-12 14:02:18
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/Core.h
 */
#pragma once
#include <string>
#include "mesh.h"
#include "SurfaceMesh.h"
namespace MeshAC{
    /* Main Functions */
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

    void generateMeshFromThreeLayer(
        const std::string &input, 
        const std::string &output, 
        double size
    );

    void adaptiveRefineMesh(
        const std::string &input,
        const std::string &output
    );


    void adaptiveRefineMeshForEdgeDislocation(
        const std::string &input,
        const std::string &output
    );









    

    /* Assistant Functions */

    void generateMeshForAtomisticRegion(
        const std::vector<Vector3D> &points,
        Mesh &resultingMesh,
        SurfaceMesh &interfaceMesh,
        bool keepConvexHull = false
    );

    void generateMeshForAtomisticRegionInThreeLayerConfiguration(
        const std::vector<Vector3D> &points,
        Mesh &resultingMesh,
        SurfaceMesh &interfaceMesh,
        Vector3D &lowerBound,
        Vector3D &upperBound,
        double &innerSize
    );

    void generateMeshWholeDomain(
        const std::vector<Vector3D> &atomPoints,
        const std::vector<Vector3D> &continnumPoints,
        Mesh &resultingMesh,
        bool keepConvexHull = false
    );

    void generateMeshForContinuumRegion(
        SurfaceMesh &boundingBoxSurfaceMesh,
        SurfaceMesh &interfaceSurfaceMesh,
        Mesh &resultingMesh,
        bool hasHoles=true
    );

    void generateMeshForBlendRegion(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        Mesh &resultingMesh
    );

    void generateMeshForBlendRegionWithEdgeDislocation(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        Mesh &resultingMesh
    );
    
    void generateTopBottomSurfaceMeshForBlendRegion(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        SurfaceMesh &resultingMesh
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
        std::vector<std::array<double, 3>> &topEdgeNodes,
        std::vector<std::array<double, 3>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges,
        std::vector<std::array<int, 2>> &bottomEdges,
        double eps=1e-13
    );

    void extractTopBottomTrianglesInSurfaceMesh(
        SurfaceMesh &surfaceMesh,
        double top,
        double bottom,
        std::vector<std::array<double, 3>> &topEdgeNodes,
        std::vector<std::array<double, 3>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges,
        std::vector<std::array<int, 2>> &bottomEdges,
        double eps=1e-13
    );

    void generateBoundingBoxSurfaceMeshWithXYHoles(
        const Vector3D &lowerBound,
        const Vector3D &upperBound,
        double size,
        std::vector<std::array<double, 3>> &topEdgeNodes,
        std::vector<std::array<double, 3>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges, 
        std::vector<std::array<int, 2>> &bottomEdges,
        SurfaceMesh &resultingSurfaceMesh,
        bool hasHoles
    );

    void generateBoundingBoxSurfaceMeshInThreeLayerConfiguration(
        const Vector3D &OuterLowerBound,
        const Vector3D &OuterUpperBound,
        const Vector3D &InnerLowerBound,
        const Vector3D &InnerUpperBound,
        double outerSize,
        double innerSize,
        SurfaceMesh &resultingSurfaceMesh
    );


    void removeIntersectionElements(
        Mesh &originalMesh,
        Mesh &expandedAtomisticMesh,
        SurfaceMesh &interfaceSurfaceMesh,
        int outerLayers=0
    );
    
    void extractBorderingSurfaceMesh(
        std::vector<Tetrahedron *>&tets, 
        SurfaceMesh &aSurface
    );



}