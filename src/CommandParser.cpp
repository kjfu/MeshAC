/*
 * @Author: Kejie Fu
 * @Date: 2023-04-10 08:43:17
 * @LastEditTime: 2023-05-23 22:06:55
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/CommandParser.cpp
 */
#include "CommandParser.h"
#include <iostream>
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
            else if (str=="-rr"){
                res.type = FT_ADAPTIVE_REFINEMENT_FOR_EDGE_DISLOCATION;
            }
            else if (str=="-hd"){
                res.type = FT_GENERATION_FROM_EDGE_DISLOCATION_POINTS;
            }
            // else if (str=="-dev"){
            //     res.type = FT_GENERATION_DEV;
            // }
            else if (str=="-i"){
                i++;
                res.input = std::string(argv[i]);
            }
            else if (str=="-o"){
                i++;
                res.output = std::string(argv[i]);
            }
            else if (str=="-h" || str=="--help"){
                std::cout << "Usage: ./MeshAC [OPTIONS] input [output]                                                                " << std::endl;
                std::cout                                                                                                               << std::endl;
                std::cout << "Postionals:                                                                                             " << std::endl;
                std::cout << "    input  TEXT REQUIRED    Input a file of initial mesh(.mesh) for mesh generation or middle files     " << std::endl;
                std::cout << "                            (.mesh, .remesh and .value) for mesh adaptation.(string, required)          " << std::endl;
                std::cout << "    output TEXT             Output a file of resulting mesh (.mesh) or a file of interpolation solutions" << std::endl;
                std::cout << "                            (.value).                                                                   " << std::endl;
                std::cout                                                                                                               << std::endl;
                std::cout << "Options:                                                                                                " << std::endl;
                std::cout << "    -h			          Print this help message and exit.                                           " << std::endl;
                std::cout << "    -i TEXT REQUIRED        Input a file of initial mesh(.mesh) for mesh generation or middle files     " << std::endl;
                std::cout << "                            (.mesh, .remesh and .value) for mesh adaptation.(string, required)          " << std::endl;
                std::cout << "    -o TEXT                 Output a file of resulting mesh (.mesh) or a file of interpolation solutions" << std::endl;
                std::cout << "                            (.value).                                                                   " << std::endl;
                std::cout << "    -s  FLOAT               Input the max sizing value for mesh generation.                             " << std::endl;
                std::cout << "    -hd                     Generate a mesh with edge dislocation.                                      " << std::endl;
                std::cout << "    -r                      Refine an existing mesh adaptively.                                         " << std::endl;
                std::cout << "    -rr                     Refine an existing mesh with edge dislocation adaptively.                   " << std::endl;

                exit(0);
            }
        }
        return 0;
    }
}