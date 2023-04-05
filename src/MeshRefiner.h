/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:42:11
 * @LastEditTime: 2023-03-28 14:11:44
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/MeshRefiner.h
 */
#pragma once
#include "mesh.h"
#include <memory>
class MeshRefiner{
public:
    Mesh* mesh;
    MeshRefiner(Mesh *m=nullptr): mesh(m){};

    void refine(std::vector<Tetrahedron*> &elements);

    void refine(int firstIndex, std::vector<int> elementIndices, int subdomain, int nodeLabel);
};