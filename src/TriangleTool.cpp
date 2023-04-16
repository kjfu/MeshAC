#include "TriangleTool.h"
#include <string.h>
namespace MeshAC{
    void generateRectangleEdges(
        std::array<double, 2> maxPos, 
        std::array<double,2> minPos, 
        double size, 
        std::vector<std::array<double,2>> &edgeNodes, std::vector<std::array<int, 2>> &edges){

        int basebase = edgeNodes.size();
        int base = edgeNodes.size();
        double dx = maxPos[0]-minPos[0];
        double dy = maxPos[1]-minPos[1];
        int numSegmentsX = dx/size;
        int numSegmentsY = dy/size;
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ minPos[0] + dx/double(numSegmentsX)*i, minPos[1]});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ maxPos[0], minPos[1]+dy/double(numSegmentsY)*i});
            edges.push_back({base+i, base+i+1});
        }	
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ maxPos[0] - dx/double(numSegmentsX)*i, maxPos[1]});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ minPos[0], maxPos[1]-dy/double(numSegmentsY)*i});
            edges.push_back({base+i, base+i+1});
        }
        edges.back()[1] = basebase;
    }

    void generateTRIANGULATEIOWithEdges(
        std::vector<std::array<double,2>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,2>> holes, 
        double maxAreaSize,  
        triangulateio &triOut
    ){
        struct triangulateio triIn;
        setNullToTRIANGULATEIO(triIn);
        setNullToTRIANGULATEIO(triOut);
        
        triIn.numberofpoints = planeNodes.size();
        triIn.numberofpointattributes= 0;
        triIn.pointlist = new double[triIn.numberofpoints*2];
        for(int i=0; i<planeNodes.size(); i++){
            triIn.pointlist[2*i] = planeNodes[i][0];
            triIn.pointlist[2*i+1] = planeNodes[i][1];
        }

        triIn.numberofsegments = edges.size();
        triIn.segmentlist = new int[triIn.numberofsegments*2];
        for(int i=0; i<edges.size(); i++){
            triIn.segmentlist[2*i] = edges[i][0]+1;
            triIn.segmentlist[2*i+1] = edges[i][1]+1;
        }

        triIn.numberofholes=0;
        triIn.numberofregions = 0;
        if(!holes.empty()){
            triIn.numberofholes = holes.size();
            triIn.holelist = new double[triIn.numberofholes*2];
            for(int i=0; i<holes.size(); i++){
                triIn.holelist[i*2] = holes[i][0];
                triIn.holelist[i*2+1] = holes[i][1];
            }
        }


        std::string str= "pqYQD";
        if(maxAreaSize>1e-10){
            str = str + "a" + std::to_string(maxAreaSize);
        }

        char cmd[256];
        strcpy(cmd, str.c_str());
        triangulate(cmd, &triIn, &triOut, nullptr);
        deleteTRIANGULATEIOAllocatedArrays(triIn);
    }

    void deleteTRIANGULATEIOAllocatedArrays(triangulateio &io){
        delete [] io.edgelist;
        delete [] io.edgemarkerlist;
        //delete [] io.holelist;
        delete [] io.neighborlist;
        delete [] io.normlist;
        delete [] io.pointattributelist;
        delete [] io.pointmarkerlist;
        delete [] io.pointlist;
        delete [] io.regionlist;
        delete [] io.segmentlist;
        delete [] io.segmentmarkerlist;
        delete [] io.trianglearealist;
        delete [] io.triangleattributelist;
        delete [] io.trianglelist;	
    }

    void setNullToTRIANGULATEIO(triangulateio &io){
        io.edgelist = nullptr;
        io.edgemarkerlist = nullptr;
        io.holelist=nullptr;
        io.neighborlist=nullptr;
        io.normlist= nullptr;
        io.pointattributelist=nullptr;
        io.pointmarkerlist=nullptr;
        io.pointlist=nullptr;
        io.regionlist=nullptr;
        io.segmentlist=nullptr;
        io.segmentmarkerlist=nullptr;
        io.trianglearealist=nullptr;
        io.triangleattributelist=nullptr;
        io.trianglelist=nullptr;
    }
}