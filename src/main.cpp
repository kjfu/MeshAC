#include "tetgen.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <vector>
#include <array>
#include "mesher3d_io.h"
#include "mesher3d_core.h"
#include <iostream>
#include <set>
#include "mesh.h"
#include <fstream>
#include <ctime>
#include <unordered_map>
using namespace MeshAC;
int main(int argc, char *argv[]){

    tetgenio in, out;
	Mesh goalMesh;
	Mesh backgroundMesh;

    char inFilePath[1024];
	int choice = -1;

    std::string str;
    REAL size = 0;
	double hmax = -1;
	double hmin = -1;

	std::string refineFileHeadIn;
    std::string refineFileHeadOut = "out3d";

	std::string fileIn;
    std::string fileOut = "out3d.mesh";
	bool quiet = false;
	for(int i=1; i<argc; i++){
        str = std::string(argv[i]);

        if (str == "-s"){
            i++;
            size = atof(argv[i]);
        }
		else if(str=="-q"){
			quiet = true;
		}
        else if (str == "-o"){
            i++;
            fileOut = std::string(argv[i]);

        }
		else if (str == "-r"){
			i++;
			refineFileHeadIn = argv[i];
			refineFileHeadOut = std::string(argv[i])+"_out";
			choice = 3;
		}
		else if(str == "-rr"){
			i++;
			refineFileHeadIn = argv[i];
			refineFileHeadOut = std::string(argv[i])+"_out";
			choice = 4;			
		}
		else if(str=="-hmin"){
            i++;
            hmin = atof(argv[i]);
		}
		else if(str=="-hmax"){
            i++;
            hmax = atof(argv[i]);
		}
		else if(str == "-hd"){
			i++;
			fileIn = argv[i];
			choice = 5;
		}
		else if(str == "-ch"){
			i++;
			fileIn = argv[i];
			choice = 6;
		}
		else if(str == "-hh"){
			i++;
			fileIn = argv[i];
			choice = 7;
		}
		else if(str == "-hd2"){
			i++;
			fileIn = argv[i];
			choice = 8;
		}
        else if (str.length()>2){
            
			if (str.find(".poly") != std::string::npos){
				size_t cut = str.find(".poly");
				strcpy(inFilePath, str.substr(0, cut).c_str());
				in.load_poly(inFilePath);
				choice = choice>0?choice:1;
			}
			else if (str.find(".mesh") != std::string::npos){
				// size_t cut = str.find(".mesh");
				// strcpy(inFilePath, str.substr(0, cut).c_str());
				// in.load_medit(inFilePath, 1);
				if (choice==3||choice==4){
					continue;
				}				
				fileIn = argv[i];
				choice = choice>0?choice:2;

			}
        }

    }

	//TEST
	//  choice = 3;
	//  quiet = true;
	// size = 5;
	choice = 9;
	fileIn = "/home/kjfu/research/Mesher3DForSJTU/examples/quality_test/CPC-disloc/step4.mesh";
	fileOut = fileIn + ".csv";
	// fileOut = "/home/kjfu/research/Mesher3DForSJTU/examples/paper_test/out3d.mesh";
	// fileIn = "C:\\Users\\Fu Kejie\\Desktop\\workspace\\Mesher3DForSJTU\\examples\\periodic\\dislocation.mesh";
	// refineFileHeadIn = "/home/kjfu/research/Mesher3DForSJTU/examples/refine_case/20210430/out3d";
	// refineFileHeadOut = "/home/kjfu/research/Mesher3DForSJTU/examples/refine_case/test3d_out";
	//TEST
    if (choice == -1){
		fprintf(stderr, "No input file!\n");
		exit(1);
    }
	if (choice == 1){
			constrainedTetrahedralization(&in, &out, size, quiet);    
			saveAsMESH(&out, fileOut);
	}
	else if (choice == 2){
		std::vector<int> tetMarkers;

		delaunayTetrahedralization(fileIn, fileOut, size, quiet);

	}
	else if (choice == 3){
		refineMeshV2(refineFileHeadIn, refineFileHeadOut, quiet);
	}
	else if (choice == 4){
		refineMeshV3(refineFileHeadIn, refineFileHeadOut, hmax, hmin);
	}
	else if(choice == 5){
		generateZHandleMeshV2(fileIn, fileOut, size, quiet);
	}
	else if(choice == 6){
		generateConvexHull(fileIn, fileOut);
	}
	else if (choice == 7){
		delaunayTetrahedralizationWithHoles(fileIn, fileOut, size, quiet);
	}
	else if(choice == 8){
		generateZHandleMeshV3(fileIn, fileOut, size, quiet);
	}
	else if(choice == 9){
		analysisQuality(fileIn, fileOut);
	}
    return 0;
}
