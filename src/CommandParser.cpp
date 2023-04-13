/*
 * @Author: Kejie Fu
 * @Date: 2023-04-10 08:43:17
 * @LastEditTime: 2023-04-13 21:12:25
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/CommandParser.cpp
 */
#include "CommandParser.h"

namespace MeshAC{
    int parseCommandLine(int argc, char* argv[], CommandInfo &res){
        for(int i=1; i<argc; i++){
            std::string str = std::string(argv[i]);
            if (str=="-s"){
                i++;
                res.size = atof(argv[i]);
            }
            else if (str=="-r"){
                res.type = FT_ADAPTIVE_REFINEMENT;
            }
            else if (str=="-hd"){
                res.type = FT_GENERATION_WITH_EDGE_DISLOCATION_POINTS;
            }
            else if (str=="-i"){
                i++;
                res.input = std::string(argv[i]);
            }
            else if (str=="-o"){
                i++;
                res.output = std::string(argv[i]);
            }
            else if (str=="-h" || str=="--help"){
                exit(0);
            }
        }
        return 1;
    }
}