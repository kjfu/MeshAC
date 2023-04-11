/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 20:08:02
 * @LastEditTime: 2023-04-11 14:48:56
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/GroupEntity.h
 */
#pragma once
#include <vector>
#include "SubEntity.h"
namespace MeshAC{
class VolumeShell{
public:
    bool closed=true;
    std::vector<SubEdge> tetrahedronEdges;
    std::vector<Node*> nodes;
    VolumeShell(){};

};

class VolumeBall{
public:
    bool closed = false;
    std::vector<SubTriangle> tetrahedronTriangles;
    VolumeBall(){};

};

class VolumeUmbrella{
public:
    std::vector<SubTriangle> tetrahedronTriangles;
    Node *endNode;
    VolumeUmbrella(){};
};
}