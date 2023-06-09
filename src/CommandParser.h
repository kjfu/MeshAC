/*
 * @Author: Kejie Fu
 * @Date: 2023-04-10 08:40:35
 * @LastEditTime: 2023-06-07 11:47:28
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/CommandParser.h
 */
#pragma once
#include <string>
namespace MeshAC{
    enum FUNCTION_TYPE{
        FT_GENERATION_FROM_POINTS=0,
        FT_GENERATION_FROM_EDGE_DISLOCATION_POINTS,
        FT_ADAPTIVE_REFINEMENT,
        FT_ADAPTIVE_REFINEMENT_FOR_EDGE_DISLOCATION,
        FT_GENERATION_THREE_LAYER
    };
    struct CommandInfo{
        FUNCTION_TYPE type = FT_GENERATION_FROM_POINTS;
        std::string input;
        std::string output;
        double size=-1;

    };

    /**
     * @brief 
     * 
     * @param argc 
     * @param argv 
     * @param res 
     * @return int 
     */
    int parseCommandLine(int argc, char* argv[], CommandInfo &res);

}
