/*
 * @Author: Kejie Fu
 * @Date: 2023-04-10 09:38:57
 * @LastEditTime: 2023-04-14 14:26:30
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/TetgenTool.h
 */
#pragma once
#include "tetgen.h"
#include "Element.h"
#include "Vector3D.h"
#include "mesh.h"
#include "SurfaceMesh.h"
#include <vector>
namespace MeshAC
{
    

    void transportNodesToTETGENIO(
        const std::vector<Node *> &sNodes, 
        tetgenio &out
    );

    void transportPointsToTETGENIO(
        const std::vector<Vector3D> &points, 
        int label, 
        tetgenio &out
    );

    void transportSurfaceMeshToTETGENIO(
        SurfaceMesh &surfaceMesh, 
        std::vector<Vector3D> &holeCenters,
        tetgenio &out
    );

    void transportTETGENIOToMesh(
        tetgenio &in, 
        Mesh &out, 
        int withLabel=false
    );

} // namespace MeshAC
