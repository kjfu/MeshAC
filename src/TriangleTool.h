/*
 * @Author: Kejie Fu
 * @Date: 2023-04-14 09:56:34
 * @LastEditTime: 2023-04-16 17:16:27
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/TriangleTool.h
 */
#pragma once
#include <array>
#include <vector>
#include "triangle.h"
namespace MeshAC{

    /**
     * @brief 
     * 
     * @param maxPos 
     * @param minPos 
     * @param size 
     * @param edgeNodes 
     * @param edges 
     */
    void generateRectangleEdges(
        std::array<double, 2> maxPos, 
        std::array<double,2> minPos, 
        double size, 
        std::vector<std::array<double,2>> &edgeNodes, std::vector<std::array<int, 2>> &edges
    );
    
    /**
     * @brief 
     * 
     * @param planeNodes 
     * @param edges 
     * @param holes 
     * @param maxAreaSize 
     * @param triOut 
     */
    void generateTRIANGULATEIOWithEdges(
        std::vector<std::array<double,2>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,2>> holes, 
        double maxAreaSize,  
        triangulateio &triOut
    );

    /**
     * @brief 
     * 
     * @param io 
     */
    void deleteTRIANGULATEIOAllocatedArrays(triangulateio &io);

    /**
     * @brief 
     * 
     * @param io 
     */
    void setNullToTRIANGULATEIO(triangulateio &io);

}