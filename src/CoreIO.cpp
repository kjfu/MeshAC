#include "CoreIO.h"
#include <fstream>
#include <sstream>
namespace MeshAC{
    void loadREMESH(
        std::vector<int> &elements, 
        std::vector<Vector3D> &points, 
        const std::string &filePath){
        std::ifstream inFile(filePath);
        if (inFile.is_open()){
            while(inFile){
                std::string line;
                std::string keyString;
                std::getline(inFile, line);
                std::stringstream lineStream(line);
                lineStream >> keyString;
                if (keyString == "Append_points"){
                    std::getline(inFile, line);
                    lineStream.clear();
                    lineStream.str(line);
                    int numPoints;
                    lineStream >> numPoints;
                    for (int i=0; i<numPoints; i++){
                        std::getline(inFile, line);
                        std::stringstream subLineStream(line);
                        double x, y, z;
                        subLineStream >> x >> y >> z;
                        points.emplace_back(x, y, z);
                    }
                }
                else if(keyString == "Refine_elements"){
                    std::getline(inFile, line);
                    lineStream.clear();
                    lineStream.str(line);
                    int numElements;
                    lineStream >> numElements;
                    for (int i=0; i<numElements; i++){
                        std::getline(inFile, line);
                        std::stringstream subLineStream(line);
                        int id;
                        subLineStream >> id;
                        elements.push_back(id);
                    }
                }

            }
        }
    
    }
}