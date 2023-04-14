/*
 * @Author: Kejie Fu
 * @Date: 2023-04-14 17:10:51
 * @LastEditTime: 2023-04-14 17:11:38
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/CoreIO.h
 */
#pragma once
#include<vector>
#include<Vector3D.h>
#include<string>
namespace MeshAC {
    void loadREMESH(
        std::vector<int> &elements, 
        std::vector<Vector3D> &points, 
        const std::string &filePath
    );
}

