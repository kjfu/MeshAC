#pragma once
#include "Vector3D.h"
#include <unordered_map>
#include <vector>

namespace MeshAC{
    class PointSet{
    public:
        std::unordered_map<int, std::vector<Vector3D> > subsets;
        PointSet(){}

        int loadPointSet(const std::string &filePath);

        int calculateBoundingBox(Vector3D &lowerBound, Vector3D &upperBound);
        int calculateSubsetBoundingBox(int label, Vector3D &lowerBound, Vector3D &upperBound);
        
    };
}