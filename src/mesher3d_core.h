/*
 * @Author: Kejie Fu
 * @Date: 2021-09-07 21:09:25
 * @LastEditTime: 2023-03-31 15:20:40
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/mesher3d_core.h
 */
#pragma once
#include "tetgen.h"
#include <vector>
#include <string>
#include "triangle.h"
#include "surfaceMesh.h"

void analysisQuality(const std::string &fileIn, const std::string &fileOut);

void generateConvexHull(const std::string &fileIn, const std::string &fileOut);
void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers, bool beQuiet=false);//useless
void delaunayTetrahedralization(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);
void delaunayTetrahedralizationWithHoles(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);

void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size, bool beQuiet=false);


void refineMesh(const std::string &fileInHead, const std::string &fileOutHead, bool beQuiet=false);
void refineMeshV2(const std::string &fileInHead, const std::string &fileOutHead, bool beQuiet=false);
void refineMeshV3(const std::string &fileInHead, const std::string &fileOutHead, double hmax, double hmin);

void generateZHandleMesh(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);
void parseZHandle(SurfaceMesh &zHandleSurface, double top, double bottom, 
	std::vector<std::array<double, 2>> &topEdgeNodes,	
    std::vector<std::array<double, 2>> &bottomEdgeNodes,
    std::vector<std::array<int, 2>> &topEdges, 
    std::vector<std::array<int, 2>> &bottomEdges);

void parseZHandleV2(SurfaceMesh &zHandleSurface, Vector3D xyzmax, Vector3D xyzmin, Vector3D oxyzmax, Vector3D oxyzmin, double size, Mesh &outMesh);
void parseZHandleV3(SurfaceMesh &zHandleSurface, Vector3D xyzmax, Vector3D xyzmin, Vector3D oxyzmax, Vector3D oxyzmin, double size, Mesh &outMesh);
void generateZHandleMeshV2(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);
void generateZHandleMeshV3(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);

void generatePeriodicBoundaryConditionMesh(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet=false);//TODO



void setNullToTRIANGULATEIO(triangulateio &io);
void deleteTRIANGULATEIOAllocatedArrays(triangulateio &io);
void generateRectangle(std::array<double, 2> maxPos, std::array<double,2> minPos, double size, std::vector<std::array<double,2>> &edgeNodes, std::vector<std::array<int, 2>> &edges);
void generateMeshInPlaneWithEdges(std::vector<std::array<double,2>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,2>> holes, double maxAreaSize, triangulateio &triOut);

//Useless!
void generateConvaxHullFromPointsInPlane(std::vector<std::array<double, 2>> &planeNodes, std::array<double, 2> &oxymax, std::array<double, 2> &oxymin, std::vector<std::array<double,2>> &finalNodes, std::vector<std::array<int, 2>> &finalEdges);

void generateConvaxHullFromPointsIn3D(tetgenio &tet, Vector3D &oxyzmax, Vector3D &oxyzmin, Mesh &goalMesh, SurfaceMesh &goalSurface, bool withPointLabel=false);
void generateConvaxHullFromPointsIn3D(tetgenio &tet, Mesh &goalMesh, SurfaceMesh &goalSurface);
void generateConvaxHullFromPointsIn3DRemoveHoles(tetgenio &tet, Mesh &goalMesh, SurfaceMesh &goalSurface);
void resetPoints(tetgenio &tet, Vector3D pMax, Vector3D pMin, std::vector<int> &indexOf1);

void extractBorder(std::vector<Tetrahedron *>&tets, SurfaceMesh &aSurface);
void generateBoundingBoxTETGENIO(tetgenio &tetIn, Vector3D pMax, Vector3D pMin, double size, tetgenio &tetOut);
void generateBoundingBoxTETGENIO(Vector3D xyzmax, Vector3D xyzmin, double size, SurfaceMesh &aSurface);