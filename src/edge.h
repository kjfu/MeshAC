#pragma once
#include <array>
#include <vector>
#include <algorithm>
#include "node.h"
class TriangleElement;
class TetrahedronElement;

class Edge
{
public:
    TriangleElement* tri = nullptr;
    TetrahedronElement* tet = nullptr;
    void* tmp;
    std::array<int, 2> orderedNodeIndices;
    std::array<Node*, 2> sNodes;
    int localIndexOfTriangle;
    int localIndexOfTetrahedron;
    
    Edge() = default;

    Edge(int i0, int i1){
        if (i0<i1){
            orderedNodeIndices[0] = i0;
            orderedNodeIndices[1] = i1;
        }
        else{
            orderedNodeIndices[0] = i1;
            orderedNodeIndices[1] = i0;
        }
    }

    Edge(Node* n0, Node* n1){
        sNodes[0] = n0;
        sNodes[1] = n1;
    }

    double length(){
        return (sNodes[0]->pos - sNodes[1]->pos).norm();
    }

    void sortNodeIndices();

};