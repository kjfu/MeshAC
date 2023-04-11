/*
 * @Author: Kejie Fu
 * @Date: 2022-03-26 22:06:24
 * @LastEditTime: 2023-04-10 08:57:29
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/mesher3d_io.h
 */
#pragma once

#include "tetgen.h"
#include "triangle.h"

#include <string>
#include <vector>
#include <array>
#include "Vector3D.h"


namespace MeshAC{

void loadMesh(tetgenio *in, std::string filePath);


void loadREMESH(std::vector<int> &elements, std::vector<std::array<double,3>> &points, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers);



void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin, std::vector<int> &indexOf1);
void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin);//useless
void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min);

}


