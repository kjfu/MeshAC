/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:42:11
 * @LastEditTime: 2023-04-10 08:52:07
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/MeshRefiner.h
 */
#pragma once
#include "mesh.h"
#include <memory>
namespace MeshAC{
class MeshRefiner{
public:
    Mesh* mesh;
    MeshRefiner(Mesh *m=nullptr): mesh(m){};

    void refine(std::vector<Tetrahedron*> &elements);

    void refine(int firstIndex, std::vector<int> elementIndices, int subdomain, int nodeLabel);
};
}