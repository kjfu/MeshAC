#pragma once
#include <vector>
#include "vector3d.h"
struct Node
{

    Vector3D pos;
    int index;
    int label=0;
    bool fixed = false;
    double sizing;
    std::vector<double> scalarValues;
    std::vector<Vector3D> vectorValues;
    int edit = 0;

    void *tempDate = nullptr;
    
    Node(){

    }

    Node(double x, double y, double z):pos(x,y,z){
    }

    Node(double *xyz): pos(xyz){
    }
    Node(Vector3D &vec):pos(vec){
    }

};