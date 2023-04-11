/*
 * @Author: Kejie Fu
 * @Date: 2023-04-06 01:09:29
 * @LastEditTime: 2023-04-11 14:20:11
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/SurfaceMesh.h
 */
#pragma once


#include <vector>
#include <string>
#include <unordered_map>
#include "tetgen.h"
#include "Element.h"
#include "triangle.h"
namespace MeshAC{
enum PROJECTION_TYPE{
    XY_PLANE = 0,
    YZ_PLANE = 1,
    ZX_PLANE = 2,
};
class SurfaceMesh{
public:
    std::vector<Node *> nodes;
    std::vector<Triangle *> triangles;
    std::unordered_map<int, std::string> zoneNames;
    SurfaceMesh(){}
    ~SurfaceMesh(){

        for(int i=0; i<nodes.size(); i++){
            delete nodes[i];
        }
        for(int i=0; i<triangles.size(); i++){
            delete triangles[i];
        }     
    }


    void cloneSurfaceMesh(SurfaceMesh &another);

    void addSubSurfaceMesh(SurfaceMesh &another);//Ensure no intersect elements with current mesh


    void mergeSurfaceMesh(SurfaceMesh &another, double tolerance = std::numeric_limits<double>::epsilon());

    void rebuildTriangleAdjacency();


    void deleteTriangles(std::vector<Triangle *> &toDelTriangles);

    
    void projectTRIANGULATEIO(triangulateio &in, PROJECTION_TYPE projectionType, double offset);

    //
    void loadMESH(const std::string &filePath);
    void loadVTK(const std::string &filePath);
    void loadPLY(const std::string &filePath);
    void loadTETGENIO(tetgenio &in);
    void exportPLY(const std::string &filePath);
    void exportVTK(const std::string &filePath);
    void exportMESH(const std::string &filePath);
    void exportSU2(const std::string &filePath);
    void exportTETGENIO(tetgenio &out);
    void rebuildIndices();

    double maxSizing = std::numeric_limits<double>::min();
    double minSizing = std::numeric_limits<double>::max();
    void estimateSizing();
};

}