/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 20:08:02
 * @LastEditTime: 2023-03-28 12:11:12
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/GroupEntities.h
 */
#pragma once
#include <vector>
#include "SubEntities.h"
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