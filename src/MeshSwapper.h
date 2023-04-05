/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:35:10
 * @LastEditTime: 2023-03-28 13:06:56
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/MeshSwapper.h
 */
#pragma once
#include "mesh.h"
#include "SubEntities.h"
#include "GroupEntities.h"
#include <memory>
class MeshSwapper{
public:
    Mesh *mesh;

    MeshSwapper(Mesh* m=nullptr):mesh(m){}

    void SwapOptimization(int SubdomainNumber,  double minQuality);

    bool checkEdgeSwap(Tetrahedron *tet, int iLocal, VolumeShell &aShell, VolumeUmbrella &anUmbrella, std::vector<double> &qualities);

    void EdgeSwap(VolumeShell &aShell,VolumeUmbrella &anUmbrella, std::vector<double> &qualities);


};