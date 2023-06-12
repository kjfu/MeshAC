#include "TriangleTool.h"
#include <string.h>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <functional>
namespace MeshAC{

    void generateRectangleEdgesWithDifferentSize(
        std::array<double, 2> maxPos, 
        std::array<double,2> minPos, 
        double leftsize, 
        double rightsize,
        double topSize,
        double bottomSize,
        std::vector<std::array<double,3>> &edgeNodes, std::vector<std::array<int, 2>> &edges
    ){
        int basebase = edgeNodes.size();
        int base = edgeNodes.size();
        double dx = maxPos[0]-minPos[0];
        double dy = maxPos[1]-minPos[1];
        int numSegmentsX = std::fmax(1, dx/bottomSize);
        numSegmentsX = numSegmentsX%2 ? numSegmentsX: numSegmentsX+1;
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ minPos[0] + dx/double(numSegmentsX)*i, minPos[1], 0});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();
        int numSegmentsY = std::fmax(1, dy/rightsize);
        numSegmentsY = numSegmentsY%2 ? numSegmentsY: numSegmentsY+1;
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ maxPos[0], minPos[1]+dy/double(numSegmentsY)*i, 0});
            edges.push_back({base+i, base+i+1});
        }	
        base = edgeNodes.size();
        numSegmentsX = std::fmax(1, dx/topSize);
        numSegmentsX = numSegmentsX%2 ? numSegmentsX: numSegmentsX+1;
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ maxPos[0] - dx/double(numSegmentsX)*i, maxPos[1], 0});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();

        numSegmentsY = std::fmax(1, dy/leftsize);
        numSegmentsY = numSegmentsY%2 ? numSegmentsY: numSegmentsY+1;
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ minPos[0], maxPos[1]-dy/double(numSegmentsY)*i, 0});
            edges.push_back({base+i, base+i+1});
        }
        edges.back()[1] = basebase;
    }



    void generateRectangleEdges(
        std::array<double, 2> maxPos, 
        std::array<double,2> minPos, 
        double size, 
        std::vector<std::array<double,3>> &edgeNodes, std::vector<std::array<int, 2>> &edges){

        int basebase = edgeNodes.size();
        int base = edgeNodes.size();
        double dx = maxPos[0]-minPos[0];
        double dy = maxPos[1]-minPos[1];
        int numSegmentsX = std::fmax(1, dx/size);
        int numSegmentsY = std::fmax(1, dy/size);
        numSegmentsX = numSegmentsX%2 ? numSegmentsX: numSegmentsX+1;
        numSegmentsY = numSegmentsY%2 ? numSegmentsY: numSegmentsY+1;
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ minPos[0] + dx/double(numSegmentsX)*i, minPos[1], 0});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ maxPos[0], minPos[1]+dy/double(numSegmentsY)*i, 0});
            edges.push_back({base+i, base+i+1});
        }	
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsX; i++){
            edgeNodes.push_back({ maxPos[0] - dx/double(numSegmentsX)*i, maxPos[1], 0});
            edges.push_back({base+i, base+i+1});
        }
        base = edgeNodes.size();
        for(int i=0; i<numSegmentsY; i++){
            edgeNodes.push_back({ minPos[0], maxPos[1]-dy/double(numSegmentsY)*i, 0});
            edges.push_back({base+i, base+i+1});
        }
        edges.back()[1] = basebase;
    }

    void generateTRIANGULATEIOWithEdges(
        std::vector<std::array<double,3>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,3>> holes, 
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

    void findHoles( std::vector<std::array<double,3>> &holeNodes, 
    std::vector<std::array<int, 2>> &holeEdges, 
    std::vector<std::vector<std::array<double,3>>> &innerLoopNodes,
    std::vector<std::vector<std::array<int, 2>>> &innerLoopEdges,
    std::vector<std::array<double,3>> &holes){
        std::vector<int> fa(holeNodes.size());
        for(int i=0; i<fa.size(); i++){
            fa[i]=i;
        }
        std::function<int(int)> findf;
        findf=[&](int i){
            if(fa[i]==i){
                return i;
            }
            return findf(fa[i]);
        };

        auto merge=[&](int i, int j){
            fa[findf(i)]=findf(j);
        };

        for(auto edg: holeEdges){
            merge(edg[0], edg[1]);
        }
        std::unordered_map<int, std::vector<int>> groupHoles;
        for(int i=0; i<fa.size(); i++){
            groupHoles[findf(i)].push_back(i);
        }
        //TODO:
        std::vector<std::array<double, 3>> upperBounds;
        std::vector<std::array<double, 3>> lowerBounds;
        std::vector<std::vector<int>> groupNodes;
        for(auto &group: groupHoles){
            std::vector<int> &nodes = group.second;
            groupNodes.push_back(nodes);
            upperBounds.push_back({holeNodes[nodes[0]][0], holeNodes[nodes[0]][1], 0});
            lowerBounds.push_back({holeNodes[nodes[0]][0], holeNodes[nodes[0]][1], 0});
            for(int i=1; i<nodes.size(); i++){
                for(int j=0; j<2; j++){
                    upperBounds.back()[j]=fmax(upperBounds.back()[j], holeNodes[nodes[i]][j]);
                    lowerBounds.back()[j]=fmin(lowerBounds.back()[j], holeNodes[nodes[i]][j]);
                }
            }
        }
        std::vector<bool> usefulHoles(groupHoles.size(), true);
        for(int i=0; i<upperBounds.size(); i++){
            if (usefulHoles[i]==false) continue;
            for(int j=0; j<lowerBounds.size(); j++){
                if(i==j) continue;
                if(usefulHoles[j]==false) continue;
                if (upperBounds[i][0]>upperBounds[j][0] && upperBounds[i][1]>upperBounds[j][1] && lowerBounds[i][0]<lowerBounds[j][0] && lowerBounds[i][1]<lowerBounds[j][1]){
                    usefulHoles[j] = false;
                }
            }
        }
        std::unordered_map<int, int> oldNodesMapNewNodes;
        std::vector<std::array<double,3>> updateHoleNodes;
        std::vector<std::array<int, 2>> updateHoleEdges;
        int startIndex=0;
        for(int i=0; i<usefulHoles.size(); i++){
            if(usefulHoles[i]){
                holes.push_back({0.5*upperBounds[i][0]+0.5*lowerBounds[i][0], 0.5*upperBounds[i][1]+0.5*lowerBounds[i][1], 0});
                for(auto v: groupNodes[i]){
                    updateHoleNodes.push_back(holeNodes[v]);
                    oldNodesMapNewNodes[v]= startIndex;
                    startIndex++;
                }
            }
            else{
                innerLoopNodes.emplace_back();
                innerLoopEdges.emplace_back();
                std::unordered_map<int, int> tmpUpdateHoleNodes;
                int tmpIndex=0;
                for(auto v : groupNodes[i]){
                    tmpUpdateHoleNodes[v] = tmpIndex;
                    tmpIndex++;
                    innerLoopNodes.back().push_back(holeNodes[v]);

                }
                for(auto e: holeEdges){
                    if (tmpUpdateHoleNodes.count(e[0])){
                        innerLoopEdges.back().push_back({tmpUpdateHoleNodes[e[0]], tmpUpdateHoleNodes[e[1]]});
                    }
                }
            }
        }
        for(auto e: holeEdges){
            if (oldNodesMapNewNodes.count(e[0])){
                updateHoleEdges.push_back({oldNodesMapNewNodes[e[0]], oldNodesMapNewNodes[e[1]]});
            }
        }
        holeNodes=updateHoleNodes;
        holeEdges=updateHoleEdges;;
    }
}