/*
 * @Author: Kejie Fu
 * @Date: 2023-04-13 23:29:22
 * @LastEditTime: 2023-04-15 10:10:15
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/PointSet.cpp
 */
#include "PointSet.h"
#include <fstream>
#include <sstream>
#include <cmath>
namespace MeshAC {
    int PointSet::loadPointSet(const std::string &filePath){
	    std::ifstream inFile(filePath);
        if (inFile.is_open()){

            while (inFile){		
                std::string line;
                std::string keystring;
                std::getline(inFile, line);
                std::stringstream lineStream(line);
                lineStream >> keystring;
                if (keystring == "Vertices"){
                    std::getline(inFile, line);
                    lineStream.clear();
                    lineStream.str(line);				
                    int nv;
                    lineStream >> nv;
                    for(int i=0; i<nv; i++){
                        std::getline(inFile, line);
                        lineStream.clear();
                        lineStream.str(line);
                        double x, y, z;
                        int marker;
                        lineStream >> x >> y >> z >> marker;
                        subsets[marker].emplace_back(x,y,z);
                    }
                }
            }

            inFile.close();
        }
        return 0;
    }
    
    int PointSet::calculateBoundingBox(Vector3D &lowerBound, Vector3D &upperBound){
        lowerBound.initialize(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        upperBound.initialize(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        for(auto &kv: subsets){
            for(auto &v: kv.second){
                for(int i=0; i<3; i++){
                    lowerBound[i] = fmin(lowerBound[i], v[i]);
                    upperBound[i] = fmax(upperBound[i], v[i]);
                }
            }
        }
        return 0;
    }

    int PointSet::calculateSubsetBoundingBox(int label, Vector3D &lowerBound, Vector3D &upperBound){
        std::vector<Vector3D> subset = subsets[label];
        lowerBound.initialize(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        upperBound.initialize(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
        for(auto &v: subset){
            for(int i=0; i<3; i++){
                lowerBound[i] = fmin(lowerBound[i], v[i]);
                upperBound[i] = fmax(upperBound[i], v[i]);
            }
        }
        return 0;
    }

}