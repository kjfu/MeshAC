#include "mesher3d_io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
void loadMesh(tetgenio *in, std::string filePath){

	in->mesh_dim = 3;
	in->numberofpointattributes = 0;  // no point attribute.
	in->numberofpointmtrs = 0;
	in->firstnumber = 1;
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
				in->numberofpoints = nv;
				in->pointlist = new REAL[nv*3];
				in->pointmarkerlist = new int[nv];
				for(int i=0; i<nv; i++){
					std::getline(inFile, line);
					lineStream.clear();
					lineStream.str(line);
					lineStream >> in->pointlist[3*i] >> in->pointlist[3*i+1] >> in->pointlist[3*i+2] >> in->pointmarkerlist[i];
				}
			}
            else if (keystring == "Tetrahedra"){
                std::getline(inFile, line);
                std::stringstream subLineStream(line);
                int nt;
                subLineStream >> nt;
                in->numberoftetrahedra = nt;
                in->numberofcorners = 4;
                // in->numberoftetrahedronattributes = 1;
                in->tetrahedronlist = new int[4*in->numberoftetrahedra];
                // in->tetrahedronattributelist = new REAL[in->numberoftetrahedra];
                for(int i=0; i<nt; i++){
                    std::getline(inFile, line);
                    std::stringstream subSubLineStream(line);
                    subSubLineStream >> in->tetrahedronlist[4*i] >> in->tetrahedronlist[4*i+1] >> in->tetrahedronlist[4*i+2] >> in->tetrahedronlist[4*i+3];
                    // in->tetrahedronattributelist[i] = 0;
                }
            }
		}
	}
}


void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin, std::vector<int> &indexOf1){
	tetIn.mesh_dim = 3;
	tetIn.numberofpointattributes = 0;  // no point attribute.
	tetIn.numberofpointmtrs = 0;
	tetIn.firstnumber = 1;
    max.xyz.fill(std::numeric_limits<double>::min());
    min.xyz.fill(std::numeric_limits<double>::max());
    omax.xyz.fill(std::numeric_limits<double>::min());
    omin.xyz.fill(std::numeric_limits<double>::max());





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
				tetIn.numberofpoints = nv;
				tetIn.pointlist = new REAL[tetIn.numberofpoints*3];
				tetIn.pointmarkerlist = new int[tetIn.numberofpoints];
                tetIn.pointmarkerlist = new int[tetIn.numberofpoints];
                int index = 0;
				for(int i=0; i<nv; i++){
					std::getline(inFile, line);
					lineStream.clear();
					lineStream.str(line);
                    double x, y, z;
                    int marker;
                    lineStream >> x >> y >> z >> marker;                        
                    tetIn.pointlist[3*i] = x;
                    tetIn.pointlist[3*i+1] = y;
                    tetIn.pointlist[3*i+2] = z;
                    tetIn.pointmarkerlist[i] = marker;
                    if (marker!=1){
                        omax[0] = std::max(omax[0], x);
                        omax[1] = std::max(omax[1], y);
                        omax[2] = std::max(omax[2], z);
                        omin[0] = std::min(omin[0], x);
                        omin[1] = std::min(omin[1], y);
                        omin[2] = std::min(omin[2], z);                        
                    }
                    else{
                        indexOf1.push_back(index);
                        max[0] = std::max(max[0], x);
                        max[1] = std::max(max[1], y);
                        max[2] = std::max(max[2], z);
                        min[0] = std::min(min[0], x);
                        min[1] = std::min(min[1], y);
                        min[2] = std::min(min[2], z);
                    }
                    index++;
				}
			}
		}

        inFile.close();
	}

}

void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min){
	tetIn.mesh_dim = 3;
	tetIn.numberofpointattributes = 0;  // no point attribute.
	tetIn.numberofpointmtrs = 0;
	tetIn.firstnumber = 1;
    max.xyz.fill(std::numeric_limits<double>::min());
    min.xyz.fill(std::numeric_limits<double>::max());


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
				tetIn.numberofpoints = nv-8;
				tetIn.pointlist = new REAL[tetIn.numberofpoints*3];
				tetIn.pointmarkerlist = new int[tetIn.numberofpoints];
                int index=0;
				for(int i=0; i<nv; i++){
					std::getline(inFile, line);
					lineStream.clear();
					lineStream.str(line);
                    double x, y, z;
                    int marker;
                    lineStream >> x >> y >> z >> marker;                        

                    if (marker!=1){
                        tetIn.pointlist[3*index] = x;
                        tetIn.pointlist[3*index+1] = y;
                        tetIn.pointlist[3*index+2] = z;
                        index++;                       
                    }
                    else{
                        max[0] = std::max(max[0], x);
                        max[1] = std::max(max[1], y);
                        max[2] = std::max(max[2], z);
                        min[0] = std::min(min[0], x);
                        min[1] = std::min(min[1], y);
                        min[2] = std::min(min[2], z);
                    }
				}
			}
		}

        inFile.close();
	}    
}

void loadNodesWithLabel(tetgenio &tetIn, std::string filePath, Vector3D &max, Vector3D &min, Vector3D &omax, Vector3D &omin){
	tetIn.mesh_dim = 3;
	tetIn.numberofpointattributes = 0;  // no point attribute.
	tetIn.numberofpointmtrs = 0;
	tetIn.firstnumber = 1;
    max.xyz.fill(std::numeric_limits<double>::min());
    min.xyz.fill(std::numeric_limits<double>::max());
    omax.xyz.fill(std::numeric_limits<double>::min());
    omin.xyz.fill(std::numeric_limits<double>::max());

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
				tetIn.numberofpoints = nv;
				tetIn.pointlist = new REAL[tetIn.numberofpoints*3];
				tetIn.pointmarkerlist = new int[tetIn.numberofpoints];
				for(int i=0; i<nv; i++){
					std::getline(inFile, line);
					lineStream.clear();
					lineStream.str(line);
                    double x, y, z;
                    int marker;
                    lineStream >> x >> y >> z >> marker;                        
                    tetIn.pointlist[3*i] = x;
                    tetIn.pointlist[3*i+1] = y;
                    tetIn.pointlist[3*i+2] = z;
                    if (marker==0){
                        omax[0] = std::max(omax[0], x);
                        omax[1] = std::max(omax[1], y);
                        omax[2] = std::max(omax[2], z);
                        omin[0] = std::min(omin[0], x);
                        omin[1] = std::min(omin[1], y);
                        omin[2] = std::min(omin[2], z);                        
                    }
                    else{
                        max[0] = std::max(max[0], x);
                        max[1] = std::max(max[1], y);
                        max[2] = std::max(max[2], z);
                        min[0] = std::min(min[0], x);
                        min[1] = std::min(min[1], y);
                        min[2] = std::min(min[2], z);
                    }
				}
			}
		}

        inFile.close();
	}

}



void loadREMESH(std::vector<int> &elements, std::vector<std::array<double,3>> &points, std::string filePath){
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
                    std::array<double,3> pos;
                    subLineStream >> pos[0] >> pos[1] >> pos[2];
                    points.push_back(pos);
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

void saveAsMESH(tetgenio *out, std::string filePath){
    std::ofstream outfile(filePath);
    outfile << "MeshVersionFormatted 2\n";
    outfile << "Dimension\n         3\n";
    outfile << "Vertices\n";
    outfile << out->numberofpoints << "\n";
    if (out->pointmarkerlist!=nullptr){
        for (int i=0; i<out->numberofpoints; i++){
            outfile << out->pointlist[i*3]
                    << "  " << out->pointlist[i*3+1]
                    << "  " << out->pointlist[i*3+2]
                    << "  " << out->pointmarkerlist[i] << "\n";
        }
    }
    else{
        for (int i=0; i<out->numberofpoints; i++){
            outfile << out->pointlist[i*3]
                    << "  " << out->pointlist[i*3+1]
                    << "  " << out->pointlist[i*3+2]
                    << "  " << 0 << "\n";
        }        
    }



    outfile << "Tetrahedra\n";
    outfile << out->numberoftetrahedra << "\n";
    for (int i=0; i<out->numberoftetrahedra; i++){
        outfile << out->tetrahedronlist[i*4]
                << "  " << out->tetrahedronlist[i*4+1]
                << "  " << out->tetrahedronlist[i*4+2]
                << "  " << out->tetrahedronlist[i*4+3]
                << "   0\n";
    }

    outfile << "End\n";
    outfile.close();

}

void saveAsMESH(tetgenio *out, std::string filePath, std::vector<int> tetMarkers){
    std::ofstream outfile(filePath);
    outfile << "MeshVersionFormatted 2\n";
    outfile << "Dimension\n         3\n";
    outfile << "Vertices\n";
    outfile << out->numberofpoints << "\n";
    for (int i=0; i<out->numberofpoints; i++){
        outfile << out->pointlist[i*3]
                << "  " << out->pointlist[i*3+1]
                << "  " << out->pointlist[i*3+2]
                << "  " << out->pointmarkerlist[i] << "\n";
    }

    outfile << "Tetrahedra\n";
    outfile << out->numberoftetrahedra << "\n";
    for (int i=0; i<out->numberoftetrahedra; i++){
        outfile << out->tetrahedronlist[i*4]
                << "  " << out->tetrahedronlist[i*4+1]
                << "  " << out->tetrahedronlist[i*4+2]
                << "  " << out->tetrahedronlist[i*4+3]
                << "  " << tetMarkers[i] << "\n";
    }

    outfile << "End\n";
    outfile.close();	
}

