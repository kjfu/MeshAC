/*
 * @Author: Kejie Fu
 * @Date: 2023-04-14 09:56:34
 * @LastEditTime: 2023-04-14 10:08:55
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/TriangleTool.h
 */
#pragma once
#include <array>
#include <vector>
#include "triangle.h"
namespace MeshAC{
    void generateRectangleEdges(
        std::array<double, 2> maxPos, 
        std::array<double,2> minPos, 
        double size, 
        std::vector<std::array<double,2>> &edgeNodes, std::vector<std::array<int, 2>> &edges
    );
    
    void generateTRIANGULATEIOWithEdges(
        std::vector<std::array<double,2>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,2>> holes, 
        double maxAreaSize,  
        triangulateio &triOut
    );

    void deleteTRIANGULATEIOAllocatedArrays(triangulateio &io);

    void setNullToTRIANGULATEIO(triangulateio &io);

}