/*
 * @Author: Kejie Fu
 * @Date: 2023-04-14 17:10:51
 * @LastEditTime: 2023-04-16 17:19:57
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/CoreIO.h
 */
#pragma once
#include<vector>
#include<Vector3D.h>
#include<string>
namespace MeshAC {
    
    /**
     * @brief 
     * 
     * @param elements 
     * @param points 
     * @param filePath 
     */
    void loadREMESH(
        std::vector<int> &elements, 
        std::vector<Vector3D> &points, 
        const std::string &filePath
    );

    
}

