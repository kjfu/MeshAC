/*
 * @Author: Kejie Fu
 * @Date: 2023-04-13 23:16:47
 * @LastEditTime: 2023-05-24 10:05:28
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
        SurfaceMesh &interfaceMesh
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