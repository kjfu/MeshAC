/*
 * @Author: Kejie Fu
 * @Date: 2023-04-10 09:38:57
 * @LastEditTime: 2023-04-11 14:19:59
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/TetgenTools.h
 */
#pragma once
#include "tetgen.h"
#include "Element.h"
#include <vector>
namespace MeshAC
{
    
    /**
     * @brief Transport a std::vector of nodes to tetgenio object.
     * 
     * @param sNodes  A std::vector of nodes.
     * @param out A result tetgenio object.
     */
    void transportNodesToTETGENIO(const std::vector<Node *> &sNodes, tetgenio &out);


} // namespace MeshAC
