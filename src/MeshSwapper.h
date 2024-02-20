/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:35:10
 * @LastEditTime: 2023-10-12 11:26:02
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/MeshSwapper.h
 */
#pragma once
#include "mesh.h"
#include "SubEntity.h"
#include "GroupEntity.h"
#include <memory>
namespace MeshAC{
class MeshSwapper{
public:
    Mesh *mesh;

    MeshSwapper(Mesh* m=nullptr):mesh(m){}

    void SwapOptimization(int SubdomainNumber,  double minQuality);

    bool checkEdgeSwap(Tetrahedron *tet, int iLocal, VolumeShell &aShell, VolumeUmbrella &anUmbrella, std::vector<double> &qualities);

    void EdgeSwap(VolumeShell &aShell,VolumeUmbrella &anUmbrella, std::vector<double> &qualities);


    bool checkFaceSwap(Tetrahedron *tet, int iLocal, std::vector<double> &qualities);
    void FaceSwap(Tetrahedron *tet, int iLocal, std::vector<double> &qualities);


};
}