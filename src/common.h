/*
 * @Author: Kejie Fu
 * @Date: 2022-03-26 22:06:24
 * @LastEditTime: 2023-04-14 22:26:11
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/common.h
 */
#pragma once

enum POINT_LABEL{
	ATOMISTIC_POINT=0,	
	CONTINUUM_POINT,
	ATOMISTIC_BORDER_POINT,
	CONTINUUM_OUT_BORDER_POINT,
};

enum TETRAHEDRON_LABEL{
	ATOMISTIC_TET=0,
	CONTINUUM_TET
};

// VTK Element types
enum VTK_ELEMENT_TYPE {	
	IGNORED = 0,
	VTK_VERTEX = 1,
	VTK_POLY_VERTEX = 2,
	VTK_LINE = 3,
	VTK_POLY_LINE = 4,
	VTK_TRIANGLE = 5,
	VTK_TRIANGLE_STRIP = 6,
	VTK_POLYGON = 7,
	VTK_PIXEL = 8,
	VTK_QUAD = 9,
	VTK_TETRA = 10,
	VTK_VOXEL = 11,
	VTK_HEXAHEDRON = 12,
	VTK_WEDGE = 13,
	VTK_PYRAMID = 14,
	VTK_QUADRATIC_EDGE = 21,
	VTK_QUADRATIC_TRIANGLE = 22,
	VTK_QUADRATIC_QUAD = 23,
	VTK_QUADRATIC_TETRA = 24,
	VTK_QUADRATIC_HEXAHEDRON = 25,
};

static const int TriangleEdge[3][2]={{1,2},{2,0},{0,1}};

static const int TetrahedronFacet[4][3] = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};

static const int TetrahedronEdge[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

static const int TetrahedronFacetShareEdge[6][2] = {{2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}};

static const int TetrahedronTwoFacetCommonEdge[4][4]={{-1, 0, 1, 2},
													{ 0,-1, 3, 4},
													{ 1, 3,-1, 5},
													{ 2, 4, 5, -1}}; 