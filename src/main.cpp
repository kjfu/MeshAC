/*
 * @Author: Kejie Fu
 * @Date: 2023-04-06 01:09:07
 * @LastEditTime: 2023-04-15 15:48:46
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/main.cpp
 */
#include "tetgen.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <array>
#include "mesher3d_io.h"
#include "mesher3d_core.h"
#include <iostream>
#include <set>
#include "mesh.h"
#include <fstream>
#include <ctime>
#include <unordered_map>
#include "CommandParser.h"
#include "Core.h"
using namespace MeshAC;
int main(int argc, char *argv[]){
	CommandInfo info;
	parseCommandLine(argc, argv, info);
	if(info.type == FT_GENERATION_FROM_POINTS){
		generateMeshFromPoints(info.input, info.output, info.size);
	}
	else if (info.type == FT_GENERATION_FROM_EDGE_DISLOCATION_POINTS){
		generateMeshFromEdgeDislocationPoints(info.input, info.output, info.size);
	}
	else if (info.type == FT_ADAPTIVE_REFINEMENT){
		adaptiveRefineMesh(info.input, info.output);
	}
	else if (info.type == FT_ADAPTIVE_REFINEMENT_FOR_EDGE_DISLOCATION){
		adaptiveRefineMeshForEdgeDislocation(info.input, info.output);
	}
	return 0;
}
