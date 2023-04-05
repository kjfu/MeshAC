/*
 * @Author: Kejie Fu
 * @Date: 2022-01-21 01:40:25
 * @LastEditTime: 2023-03-28 16:40:08
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/SpatialSearcher.h
 */
#pragma once
#include "AABB.h"
#include <unordered_map>
#include <stack>
#include "vector3d.h"
class Mesh;
class Tetrahedron;

enum POSITION_TYPE {
    INSIDE = 0,
    ONEDGE = 1,
    ONFACE = 2,
    ONVERTEX = 3,
    OUTSIDE = 4
};

struct SearchTetrahedronResult{
    Tetrahedron *tet = nullptr;
    std::array<double, 4> weights;
    POSITION_TYPE positionType = POSITION_TYPE::OUTSIDE;
    int localIndex;
};

class SpatialSearcher {

public:
    Mesh *mesh;
    aabb::Tree aTree;
    std::unordered_map<int, Tetrahedron *> ID2TET;
    std::unordered_map<Tetrahedron*, int> TET2ID;
    int currentId = 0;
    std::stack<unsigned int> trashIds;
    unsigned int getId();
    SpatialSearcher():aTree(3,0) {}
    SpatialSearcher(Mesh *m): mesh(m), aTree(3,0){ buildAABBTree();}
    void buildAABBTree();
    void searchTetrahedronContain(Vector3D position, SearchTetrahedronResult &result);
    void removeTetrahedron(Tetrahedron *pTet);
    void insertTetrahedron(Tetrahedron *pTet);
    void updateTetrahedron(Tetrahedron *pTet);

};
