/*
 * @Author: Kejie Fu
 * @Date: 2023-03-28 01:53:43
 * @LastEditTime: 2023-03-28 16:28:13
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/tetrahedron.h
 */
#pragma once
#include "node.h"
#include "aabbox.h"
#include "sphere.h"
#include <array>
#include "SubEntities.h"
#include "common.h"
#include "quality.h"
#include <iostream>
inline double SixTimesTetrahedronVolume(Vector3D v0, Vector3D v1, Vector3D v2, Vector3D v3){
    Vector3D a = v0 - v2;
    Vector3D b = v1 - v0;
    Vector3D c = v3 - v0;
    double rst = c[0]*(a[1]*b[2]-b[1]*a[2]) - c[1]*(a[0]*b[2]-b[0]*a[2]) + c[2]*(a[0]*b[1]-b[0]*a[1]);
    return rst;
}

class Tetrahedron 
{   

public:
    std::array<Node*, 4> nodes;
    std::array<Tetrahedron*, 4> adjacentTetrahedrons;
    double quality;
    AABBox boundingBox;
    Sphere circumsphere;
    int index;
    int label=0;
    bool fixed = false;
    int edit = 0;
    Tetrahedron(){
        adjacentTetrahedrons.fill(nullptr);
    }
    Tetrahedron(Node *n0, Node *n1, Node *n2, Node *n3){
        nodes[0]=n0; nodes[1]=n1; nodes[2]=n2; nodes[3]=n3;
        adjacentTetrahedrons.fill(nullptr);
    }

    Vector3D center(){
        Vector3D pos;
        for(auto n: nodes){
            pos += n->pos;
        }
        pos /= 4.0;
        return pos;
    }

    int NodeLocalIndex(Node *node){
        int rst = -1;
        for(int i=0; i<4; i++){
            if (nodes[i]==node){
                rst = i;
                break;
            }
        }
        return rst;
    }

    void calculateQuality(){
        quality = calculateTetrahedronQualityWith4Points_ISO(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, nodes[3]->pos);
    }

    void getBoundingBox(std::vector<double> &lowerBound, std::vector<double> &upperBound){
        lowerBound.resize(3);
        upperBound.resize(3);
        for(int i = 0; i <3; i++){
            lowerBound[i] = nodes[0]->pos[i];
            upperBound[i] = nodes[0]->pos[i];
        }
        for(int n=1; n<4; n++){
            for(int i = 0; i <3; i++){
                lowerBound[i] = fmin(lowerBound[i], nodes[n]->pos[i]);
                upperBound[i] = fmax(upperBound[i], nodes[n]->pos[i]);
            }            
        }        
    }

    int getLocalIndex(SubEdge &edg){
        int rst = -1;
        for(int i=0; i<6; i++){
            SubEdge sub = getSubEdge(i);
            if (sub.sortedForms[0]==edg.sortedForms[0] 
            && sub.sortedForms[1]==edg.sortedForms[1] ){
                rst = i;
                break;
            }
        }
        return rst;        
    }

    int getLocalIndex(SubTriangle &tri){
        int rst = -1;
        for(int i=0; i<4; i++){
            SubTriangle sub = getSubTriangle(i);
            if (sub.sortedForms[0]==tri.sortedForms[0] 
            && sub.sortedForms[1]==tri.sortedForms[1] 
            && sub.sortedForms[2]==tri.sortedForms[2] ){
                rst = i;
                break;
            }
        }
        return rst;
    }

    SubTriangle getSubTriangle(int iLocal){

        SubTriangle sub(nodes[TetrahedronFacet[iLocal][0]]
                        , nodes[TetrahedronFacet[iLocal][1]]
                        , nodes[TetrahedronFacet[iLocal][2]]);
        sub.tetrahedron = this;
        sub.iLocal = iLocal;
        return sub;
    }

    SubEdge getSubEdge(int iLocal){
        SubEdge sub(nodes[TetrahedronEdge[iLocal][0]], nodes[TetrahedronEdge[iLocal][1]]);
        sub.tetrahedron = this;
        sub.iLocal = iLocal;
        return sub;
    }

    // TriangleFacet facet(int index, bool withNodes = false){
    //     static int facetIndices[4][3]=
    //     { {1, 3, 2}
    //     , {0, 2, 3}
    //     , {0, 3, 1}
    //     , {0, 1, 2}};
    //     TriangleFacet rst(nodes[facetIndices[index][0]]->index
    //                     , nodes[facetIndices[index][1]]->index
    //                     , nodes[facetIndices[index][2]]->index);
    //     if (withNodes){
    //         rst.sNodes.push_back(nodes[facetIndices[index][0]]);
    //         rst.sNodes.push_back(nodes[facetIndices[index][1]]);
    //         rst.sNodes.push_back(nodes[facetIndices[index][2]]);
    //     }
    //     rst.localIndex = index;
    //     rst.tet = this;
    //     return rst;
    // }

    Node *nodeOnFacet(int facetIndex, int nodeIndex){
        static int facetIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        return nodes[facetIndices[facetIndex][nodeIndex]];
    }

    bool removeAdjacentTetrahedron(Tetrahedron *tet){
        bool rst = false;
        for(int i=0; i<4; i++){
            if(adjacentTetrahedrons[i] == tet){
                adjacentTetrahedrons[i] = nullptr;
                rst = true;
                break;
            }
        }
        return rst;
    }
    // bool dualFacet(Node* n, TriangleFacet &goalFacet){
    //     bool rst = true;
    //     int goalIndex = -1;
    //     for(int i=0; i<4; i++){
    //         if(nodes[i]==n){
    //             goalIndex = i;
    //             break;
    //         }
    //     }
    //     if (goalIndex>=0){
    //         goalFacet = facet(goalIndex);
    //     }
    //     else{
    //         rst = false;
    //     }
    //     return rst;
    // }

    bool replaceAdjacentTetrahedron(Tetrahedron *origin, Tetrahedron *another){
        bool rst = false;
        for(int i=0; i<4; i++){
            if(adjacentTetrahedrons[i] == origin){
                adjacentTetrahedrons[i] = another;
                rst = true;
                break;
            }
        }
        return rst;        
    }


    bool contain(Vector3D pos, double eps = std::numeric_limits<double>::epsilon()){
        static int nodeIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        bool rst = true;
        for(int i=0; i<4; i++){
            if (SixTimesTetrahedronVolume(nodes[nodeIndices[i][0]]->pos, nodes[nodeIndices[i][1]]->pos, nodes[nodeIndices[i][2]]->pos, pos)<(-eps)){
                rst = false;
                break;
            }
        }
        return rst;        
    }

    bool contain(Vector3D pos, std::array<double, 4> &weights){
        static int nodeIndices[4][3]=
        { {1, 3, 2}
        , {0, 2, 3}
        , {0, 3, 1}
        , {0, 1, 2}};
        bool rst = true;
        for(int i=0; i<4; i++){
            double weight = SixTimesTetrahedronVolume(nodes[nodeIndices[i][0]]->pos, nodes[nodeIndices[i][1]]->pos, nodes[nodeIndices[i][2]]->pos, pos);
            if (weight<-std::numeric_limits<double>::epsilon()){
                rst = false;
                break;
            }
            else{
                weights[i] = weight;
            }
        }
        return rst;        
    }

    double volume(){
        return SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, nodes[3]->pos)/6.0;
    }

    bool isValid(){
        return (SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, nodes[3]->pos)>0);
    }

    void interpolateNodeValue(Vector3D pos, std::vector<double> &scalars, std::vector<Vector3D> &vectors){

        double v0 = SixTimesTetrahedronVolume(nodes[1]->pos, nodes[3]->pos, nodes[2]->pos, pos);
        double v1 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[2]->pos, nodes[3]->pos, pos);
        double v2 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[3]->pos, nodes[1]->pos, pos);
        double v3 = SixTimesTetrahedronVolume(nodes[0]->pos, nodes[1]->pos, nodes[2]->pos, pos);
        double v= v0 + v1 + v2 + v3;
        std::array<double, 4> weights{{v0/v, v1/v, v2/v, v3/v}};
        int numScalars = nodes[0]->scalarValues.size();
        int numVectors = nodes[0]->vectorValues.size();
        for(int i=0; i<numScalars; i++){
            double value =0;
            for(int j=0; j<4; j++){
                value += weights[j] * nodes[i]->scalarValues[j];
            }
            scalars.push_back(value);
        }

        for(int i=0; i<numVectors; i++){
            Vector3D vec(0.0);
            for(int j=0; j<4; j++){
                vec += weights[j] * nodes[i]->vectorValues[j];
            }
            vectors.push_back(vec);
        }

    }

    void interpolateNodeValue(std::array<double, 4> &sWeights, std::vector<double> &scalars, std::vector<Vector3D> &vectors){
        std::array<double, 4> weights = sWeights;
        double sumWeights = 0;
        for(auto w: weights){
            sumWeights += w;
        }
        for(auto &w: weights){
            w /= sumWeights;
        }

        int numScalars = nodes[0]->scalarValues.size();
        int numVectors = nodes[0]->vectorValues.size();
        for(int i=0; i<numScalars; i++){
            double value =0;
            for(int j=0; j<4; j++){
                value += weights[j] * nodes[i]->scalarValues[j];
            }
            scalars.push_back(value);
        }

        for(int i=0; i<numVectors; i++){
            Vector3D vec(0.0);
            for(int j=0; j<4; j++){
                vec += weights[j] * nodes[i]->vectorValues[j];
            }
            vectors.push_back(vec);
        }
    }

    void generateBoundingBox(){
        boundingBox.reset();
        for(auto n: nodes){
            boundingBox.insert(n->pos.XYZ());
        }
    }

    void generateCircumSphere(){
        circumsphere.initialize(nodes[0]->pos
        , nodes[1]->pos
        , nodes[2]->pos
        , nodes[3]->pos);
    }

};

