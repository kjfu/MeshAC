/*
 * @Author: Kejie Fu
 * @Date: 2023-04-06 01:09:29
 * @LastEditTime: 2023-04-16 14:47:58
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

    Triangle * addTriangle(Node* n0, Node* n1, Node *n2){
        Triangle *tri = new Triangle(n0,n1,n2);
        triangles.push_back(tri);
        return tri;
    }

    void getBoundingBox(std::vector<double> &lowerBound, std::vector<double> &upperBound){
        lowerBound.resize(3);
        upperBound.resize(3);
        for(int i = 0; i <3; i++){
            lowerBound[i] = nodes[0]->pos[i];
            upperBound[i] = nodes[0]->pos[i];
        }
        for(int n=1; n<nodes.size(); n++){
            for(int i = 0; i <3; i++){
                lowerBound[i] = fmin(lowerBound[i], nodes[n]->pos[i]);
                upperBound[i] = fmax(upperBound[i], nodes[n]->pos[i]);
            }            
        }        
    }

    void cloneSurfaceMesh(SurfaceMesh &another);

    void addSubSurfaceMesh(SurfaceMesh &another);//Ensure no intersect elements with current mesh


    void mergeSurfaceMesh(SurfaceMesh &another, double tolerance = std::numeric_limits<double>::epsilon());

    void rebuildTriangleAdjacency();


    void deleteTriangles(std::vector<Triangle *> &toDelTriangles);

    void getSubRegionCenters(std::vector<Vector3D> &positions);
    
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