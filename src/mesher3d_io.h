/*
 * @Author: Kejie Fu
 * @Date: 2022-03-26 22:06:24
 * @LastEditTime: 2022-05-06 12:44:30
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/mesher3d_io.h
 */
#ifndef _MESHER3D_IO_
#define _MESHER3D_IO_

#include "tetgen.h"
#include "triangle.h"

#include <string>
#include <vector>
#include <array>
#include "vector3d.h"




void loadMesh(tetgenio *in, std::string filePath);


void loadREMESH(std::vector<int> &elements, std::vector<std::array<double,3>> &points, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath);

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers);



void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin, std::vector<int> &indexOf1);
void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin);//useless
void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min);


#endif
