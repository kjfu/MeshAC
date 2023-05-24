/*
 * @Author: Kejie Fu
 * @Date: 2023-04-13 23:16:56
 * @LastEditTime: 2023-05-24 10:15:51
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAC/src/Core.cpp
 */
#include "Core.h"
#include "CoreIO.h"
#include "PointSet.h"
#include "TetgenTool.h"
#include "TriangleTool.h"
#include "MeshRefiner.h"
#include "common.h"
#include "aabbox.h"
#include <unordered_set>
#include <unordered_map>
namespace MeshAC{
    void generateMeshFromPoints(
        const std::string &input, 
        const std::string &output, 
        double size
    ){
        //Input
        PointSet aPointSet;
        aPointSet.loadPointSet(input);
    
        //Step 1: Generate sub-mesh for atomistic region
        Mesh atomisticMesh;
        SurfaceMesh interfaceSurfaceMesh;
        generateMeshForAtomisticRegion(aPointSet.subsets[ATOMISTIC_POINT], atomisticMesh, interfaceSurfaceMesh);

        //Step 2: Generate sub-mesh for continuum region
        Mesh continuumMesh;
        Vector3D upperBound;
        Vector3D lowerBound;
        aPointSet.calculateSubsetBoundingBox(CONTINUUM_POINT, lowerBound, upperBound);
        SurfaceMesh boundingBoxSurfaceMesh;
        Vector3D dxyz=upperBound-lowerBound;
        double tmp = fmax(dxyz[0], fmax(dxyz[1], dxyz[2]));
        size = size<1e-9 ? 0.1*tmp: size;
        generateBoundingBoxSurfaceMesh(lowerBound, upperBound, size, boundingBoxSurfaceMesh);
        // boundingBoxSurfaceMesh.exportVTK("/home/kjfu/research/MeshAC/tmp/case7_multi_holes/bd.vtk");
        // interfaceSurfaceMesh.exportVTK("/home/kjfu/research/MeshAC/tmp/case7_multi_holes/interface.vtk");

        generateMeshForContinuumRegion(boundingBoxSurfaceMesh, interfaceSurfaceMesh, continuumMesh);
        
        //Output
        atomisticMesh.mergeMesh(continuumMesh, 1e-8);
        atomisticMesh.exportMESH(output);
        atomisticMesh.exportVTK(output+".vtk");
        
    }

    void generateMeshFromEdgeDislocationPoints(
        const std::string &input, 
        const std::string &output, 
        double size
    ){
        //Input
        PointSet aPointSet;
        aPointSet.loadPointSet(input);

        //Step 1: Generate sub-mesh for atomistic region
        Mesh atomisticMesh;
        SurfaceMesh interfaceSurfaceMesh;
        generateMeshForAtomisticRegion(aPointSet.subsets[ATOMISTIC_POINT], atomisticMesh, interfaceSurfaceMesh);

        //Step 2: Generate sub-mesh for continuum region
        Mesh continuumMesh;
        Vector3D upperBound;
        Vector3D lowerBound;
        aPointSet.calculateSubsetBoundingBox(CONTINUUM_POINT, lowerBound, upperBound);

        Vector3D dxyz=upperBound-lowerBound;
        double tmp = fmax(dxyz[0], fmax(dxyz[1], dxyz[2]));
        size = size<1e-9 ? 0.1*tmp: size;

        SurfaceMesh boundingBoxSurfaceMesh;
        std::vector<std::array<double, 3>> topEdgeNodes;
        std::vector<std::array<double, 3>> bottomEdgeNodes;
        std::vector<std::array<int, 2>> topEdges;
        std::vector<std::array<int, 2>> bottomEdges;
        tetgenmesh tetmesh;
        std::vector<double> tetRadius;
        double tmpCenter[3];
        double radius;
        for(auto &tet: atomisticMesh.tetrahedrons){
            tet->edit = 0;
            if(!tetmesh.circumsphere( tet->nodes[0]->pos.data(),tet->nodes[1]->pos.data(),tet->nodes[2]->pos.data(),tet->nodes[3]->pos.data(), tmpCenter, &radius)){
                radius = std::numeric_limits<double>::max();
            }
            tetRadius.push_back(radius);
        }
        std::vector<double> tmpTetRadius = tetRadius;
        std::sort(tmpTetRadius.begin(), tmpTetRadius.end());
        double eps = tmpTetRadius[tetRadius.size()*0.2];
        removeTopBottomTrianglesInSurfaceMesh(interfaceSurfaceMesh, upperBound[2], lowerBound[2], topEdgeNodes, bottomEdgeNodes, topEdges, bottomEdges, eps);
        // interfaceSurfaceMesh.exportVTK("/home/kjfu/research/MeshAC/tmp/case6_loop_max/interface.vtk");
        bool hasHoles;
        generateBoundingBoxSurfaceMeshWithXYHoles(lowerBound, upperBound, size, topEdgeNodes, bottomEdgeNodes, topEdges, bottomEdges, boundingBoxSurfaceMesh, hasHoles);
        // boundingBoxSurfaceMesh.exportVTK("/home/kjfu/research/MeshAC/tmp/case6_loop_max/bd_holes.vtk");
        generateMeshForContinuumRegion(
            boundingBoxSurfaceMesh,
            interfaceSurfaceMesh,
            continuumMesh,
            hasHoles
        );

        //Output
        atomisticMesh.mergeMesh(continuumMesh, 1e-8);
        atomisticMesh.exportMESH(output);
        atomisticMesh.exportVTK(output+".vtk");
    }

    void adaptiveRefineMesh(
        const std::string &input,
        const std::string &output
    ){
        //Input
        Mesh backgroundMesh;
        Mesh goalMesh;
        backgroundMesh.loadMESH(input + ".mesh");
	    backgroundMesh.loadNodeValues(input + ".value");

        goalMesh.loadMESH(input + ".mesh");
        std::vector<int> refine_elements;
        std::vector<Vector3D> append_points;
	    loadREMESH(refine_elements, append_points, input+".remesh");

        //Step 1:
        MeshRefiner aRefiner(&goalMesh);
        aRefiner.refine(1, refine_elements, 1, 3);
        
        //Step 2:
        std::vector<Vector3D> points;
        for(auto &n: goalMesh.nodes){
            if(n->label==2 || n->label==0){
                points.emplace_back(n->pos);
            }
        }
        for(auto &p: append_points){
            points.emplace_back(p[0], p[1], p[2]);
        }

        Mesh atomisticMesh;
        SurfaceMesh innerInterfaceSurfaceMesh;
        SurfaceMesh outerInterfaceSurfaceMesh;
        generateMeshForAtomisticRegion(points, atomisticMesh, innerInterfaceSurfaceMesh);
        removeIntersectionElements(goalMesh, atomisticMesh, outerInterfaceSurfaceMesh);

        Mesh blendMesh;
        generateMeshForBlendRegion(outerInterfaceSurfaceMesh, innerInterfaceSurfaceMesh, blendMesh);
        

        //Output
        atomisticMesh.mergeMesh(blendMesh, 1e-10);
        goalMesh.mergeMesh(atomisticMesh, 1e-10);
        backgroundMesh.interpolateNodeValuesForAnotherMesh(goalMesh);
        goalMesh.exportNodeValues(output + ".value");
        goalMesh.exportMESH(output + ".mesh");
        goalMesh.exportVTK(output+".vtk");
    }

    void adaptiveRefineMeshForEdgeDislocation(
        const std::string &input,
        const std::string &output
    ){
        //Input
        Mesh backgroundMesh;
        Mesh goalMesh;
        backgroundMesh.loadMESH(input + ".mesh");
	    backgroundMesh.loadNodeValues(input + ".value");
        
        goalMesh.loadMESH(input + ".mesh");
        std::vector<int> refine_elements;
        std::vector<Vector3D> append_points;
	    loadREMESH(refine_elements, append_points, input+".remesh");

        //Step 1:
        MeshRefiner aRefiner(&goalMesh);
        aRefiner.refine(1, refine_elements, 1, 3);
        
        //Step 2:
        std::vector<Vector3D> points;
        for(auto &n: goalMesh.nodes){
            if(n->label==2 || n->label==0){
                points.emplace_back(n->pos);
            }
        }
        for(auto &p: append_points){
            points.emplace_back(p[0], p[1], p[2]);
        }

        Mesh atomisticMesh;
        SurfaceMesh innerInterfaceSurfaceMesh;
        SurfaceMesh outerInterfaceSurfaceMesh;
        generateMeshForAtomisticRegion(points, atomisticMesh, innerInterfaceSurfaceMesh);

        removeIntersectionElements(goalMesh, atomisticMesh, outerInterfaceSurfaceMesh);

        Mesh blendMesh;
        generateMeshForBlendRegionWithEdgeDislocation(outerInterfaceSurfaceMesh, innerInterfaceSurfaceMesh, blendMesh);
        

        //Output
        atomisticMesh.mergeMesh(blendMesh, 1e-10);
        goalMesh.mergeMesh(atomisticMesh, 1e-10);
        backgroundMesh.interpolateNodeValuesForAnotherMesh(goalMesh);
        goalMesh.exportNodeValues(output + ".value");
        goalMesh.exportMESH(output + ".mesh");
        goalMesh.exportVTK(output+".vtk");        
    }



    void generateMeshForAtomisticRegion(
        const std::vector<Vector3D> &points,
        Mesh &resultingMesh,
        SurfaceMesh &interfaceMesh
    ){
        tetgenio atomisticIn;
        tetgenio atomisticOut;
        transportPointsToTETGENIO(points, ATOMISTIC_POINT, atomisticIn);

        char cmd[] = "Q";
	    tetrahedralize(cmd, &atomisticIn, &atomisticOut);

        transportTETGENIOToMesh(atomisticOut, resultingMesh);
        resultingMesh.rebuildTetrahedronsAdjacency();
        std::unordered_map<Node*, Node*> oldNewNodes;
        auto getNode
        =[&oldNewNodes]
        (Node* n){
            Node* rst =nullptr;
            if(oldNewNodes.find(n)!=oldNewNodes.end()){
                rst = oldNewNodes[n];
            }
            else{
                rst = new Node(n->pos);
                rst->label = n->label;
                oldNewNodes[n] = rst;
            }
            return rst;
        };

        std::vector<double> tetRadius;
        tetgenmesh tetmesh;
        double tmpCenter[3];
        double radius;
        for(auto &tet: resultingMesh.tetrahedrons){
            tet->edit = 0;
            if(!tetmesh.circumsphere( tet->nodes[0]->pos.data(),tet->nodes[1]->pos.data(),tet->nodes[2]->pos.data(),tet->nodes[3]->pos.data(), tmpCenter, &radius)){
                radius = std::numeric_limits<double>::max();
            }
            tetRadius.push_back(radius);
        }
        std::vector<double> tmpTetRadius = tetRadius;
        std::sort(tmpTetRadius.begin(), tmpTetRadius.end());
        double eps = tmpTetRadius[tetRadius.size()*0.2];

        int removes=0;
        do{
            removes=0;
            for(int i = 0; i <resultingMesh.tetrahedrons.size(); i++){
                if (resultingMesh.tetrahedrons[i]->edit) continue;
                bool isBorder = false;
                for (int j = 0; j < 4; j++){
                    if (resultingMesh.tetrahedrons[i]->adjacentTetrahedrons[j]==nullptr || resultingMesh.tetrahedrons[i]->adjacentTetrahedrons[j]->edit){
                        isBorder = true;
                        break;
                    }
                }
                if (isBorder){
                    if (tetRadius[i]>2*eps){
                        resultingMesh.tetrahedrons[i]->edit = 1;
                        removes++;
                    }
                }
            }
        }while(removes);

        for(auto &e: resultingMesh.tetrahedrons){
            if(e->edit) continue;
            for(int i=0; i<4; i++){
                if (e->adjacentTetrahedrons[i]==nullptr || e->adjacentTetrahedrons[i]->edit){
                    e->nodes[(i+1)%4]->label = ATOMISTIC_BORDER_POINT;
                    e->nodes[(i+2)%4]->label = ATOMISTIC_BORDER_POINT;
                    e->nodes[(i+3)%4]->label = ATOMISTIC_BORDER_POINT;
                    interfaceMesh.addTriangle(getNode(e->nodes[(i+1)%4]), getNode(e->nodes[(i+2)%4]), getNode(e->nodes[(i+3)%4]));
                }
            }
        }

        for(auto &n: oldNewNodes){
            interfaceMesh.nodes.push_back(n.second);
        }

        for(int i=0; i<resultingMesh.tetrahedrons.size(); i++){
            if (resultingMesh.tetrahedrons[i]->edit){
                delete resultingMesh.tetrahedrons[i];
                resultingMesh.tetrahedrons[i]=resultingMesh.tetrahedrons.back();
                resultingMesh.tetrahedrons.pop_back();
                i--;
            }
	    }

        resultingMesh.rebuildIndices();
        interfaceMesh.rebuildIndices();

    }

    void generateMeshForContinuumRegion(
        SurfaceMesh &boundingBoxSurfaceMesh,
        SurfaceMesh &interfaceSurfaceMesh,
        Mesh &resultingMesh,
        bool hasHoles
    ){
        for(auto n : interfaceSurfaceMesh.nodes){
            n->label = 2;
        }

        std::vector<Vector3D> holeCenters;

        if (hasHoles) interfaceSurfaceMesh.getSubRegionCenters(holeCenters);

        interfaceSurfaceMesh.mergeSurfaceMesh(boundingBoxSurfaceMesh, 1e-8);
        interfaceSurfaceMesh.removeUnclosedMesh();
        
        int numNodesOfInterface = 0;
        for(auto n: interfaceSurfaceMesh.nodes){
            if(n->label==10) numNodesOfInterface++;
        }

        int numNodesOfSurfaceMesh = interfaceSurfaceMesh.nodes.size();
        tetgenio in;
        tetgenio out;
        transportSurfaceMeshToTETGENIO(interfaceSurfaceMesh, holeCenters, in);
        char cmd[] = "pq1.1/10YQ";
	    tetrahedralize(cmd, &in, &out);

        transportTETGENIOToMesh(out, resultingMesh);
        for(int i=0; i<interfaceSurfaceMesh.nodes.size(); i++){
            if(interfaceSurfaceMesh.nodes[i]->label==2) {
                resultingMesh.nodes[i]->label = 2;
            }
            else{
                resultingMesh.nodes[i]->label = 1;
            } 
        }
        // for(int i=numNodesOfInterface; i<numNodesOfSurfaceMesh; i++){
        //     resultingMesh.nodes[i]->label = 1;
        // }
        for( int i=interfaceSurfaceMesh.nodes.size(); i<resultingMesh.nodes.size(); i++){
            resultingMesh.nodes[i]->label = 3;
        }
        for(auto &tet: resultingMesh.tetrahedrons){
            tet->label=1;
        }
    }

    void generateMeshForBlendRegion(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        Mesh &resultingMesh
    ){
        int numNodesOfInner = innerInterfaceSurfaceMesh.nodes.size();
        int numNodesOfOuter = outerInterfaceSurfaceMesh.nodes.size();

        std::vector<Vector3D> holeCenters;

        innerInterfaceSurfaceMesh.getSubRegionCenters(holeCenters);

        innerInterfaceSurfaceMesh.addSubSurfaceMesh(outerInterfaceSurfaceMesh);

        tetgenio in;
        tetgenio out;
        transportSurfaceMeshToTETGENIO(innerInterfaceSurfaceMesh, holeCenters, in);
        char cmd[] = "pq1.1/10YQ";
	    tetrahedralize(cmd, &in, &out);

        transportTETGENIOToMesh(out, resultingMesh);
        for(int i=0; i<numNodesOfInner; i++){
            resultingMesh.nodes[i]->label = 2;
        }

        for( int i=numNodesOfInner; i<resultingMesh.nodes.size(); i++){
            resultingMesh.nodes[i]->label = 3;
        }
        for(auto &tet: resultingMesh.tetrahedrons){
            tet->label=1;
        }        
    }

    void generateMeshForBlendRegionWithEdgeDislocation(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        Mesh &resultingMesh
    ){



        std::vector<Vector3D> holeCenters;

        innerInterfaceSurfaceMesh.getSubRegionCenters(holeCenters);

        SurfaceMesh blendSurfaceMesh;
        generateTopBottomSurfaceMeshForBlendRegion(outerInterfaceSurfaceMesh, innerInterfaceSurfaceMesh, blendSurfaceMesh);
        int numNodesOfInner = innerInterfaceSurfaceMesh.nodes.size();
        innerInterfaceSurfaceMesh.mergeSurfaceMesh(blendSurfaceMesh);
        int numNodesOfInnerBlend = innerInterfaceSurfaceMesh.nodes.size();
        
        innerInterfaceSurfaceMesh.mergeSurfaceMesh(outerInterfaceSurfaceMesh);

        tetgenio in;
        tetgenio out;
        transportSurfaceMeshToTETGENIO(innerInterfaceSurfaceMesh, holeCenters, in);
        char cmd[] = "pq1.1/10YQ";
	    tetrahedralize(cmd, &in, &out);

        transportTETGENIOToMesh(out, resultingMesh);
        for(int i=0; i<numNodesOfInner; i++){
            resultingMesh.nodes[i]->label = 2;
        }

        for( int i=numNodesOfInner; i<numNodesOfInnerBlend; i++){
            resultingMesh.nodes[i]->label = 1;
        }

        for( int i=numNodesOfInnerBlend; i<resultingMesh.nodes.size(); i++){
            resultingMesh.nodes[i]->label = 3;
        }
        for(auto &tet: resultingMesh.tetrahedrons){
            tet->label=1;
        }
    }


    void generateTopBottomSurfaceMeshForBlendRegion(
        SurfaceMesh &outerInterfaceSurfaceMesh,
        SurfaceMesh &innerInterfaceSurfaceMesh,
        SurfaceMesh &resultingMesh
    ){
        std::vector<double> lowerBound;
        std::vector<double> upperBound;
        innerInterfaceSurfaceMesh.getBoundingBox(lowerBound, upperBound);

        std::vector<std::array<double, 3>> outerTopEdgeNodes;
        std::vector<std::array<double, 3>> outerBottomEdgeNodes;
        std::vector<std::array<int, 2>> outerTopEdges;
        std::vector<std::array<int, 2>> outerBottomEdges;

        removeTopBottomTrianglesInSurfaceMesh(outerInterfaceSurfaceMesh, upperBound[2], lowerBound[2], outerTopEdgeNodes, outerBottomEdgeNodes, outerTopEdges, outerBottomEdges);

        std::vector<std::array<double, 3>> innerTopEdgeNodes;
        std::vector<std::array<double, 3>> innerBottomEdgeNodes;
        std::vector<std::array<int, 2>> innerTopEdges;
        std::vector<std::array<int, 2>> innerBottomEdges;
        removeTopBottomTrianglesInSurfaceMesh(innerInterfaceSurfaceMesh, upperBound[2], lowerBound[2], innerTopEdgeNodes, innerBottomEdgeNodes, innerTopEdges, innerBottomEdges);

        std::vector<std::array<double,3>> topHoles;
        std::array<double,3> topHole={0,0,0};
        double topFactor = 1.0/ double(innerTopEdgeNodes.size());
        for(auto &n: innerTopEdgeNodes){
            topHole[0]+=topFactor*n[0];
            topHole[1]+=topFactor*n[1];
        }
        topHoles.push_back(topHole);

        int base=outerTopEdgeNodes.size();
        for(auto &p: innerTopEdgeNodes){
            outerTopEdgeNodes.push_back(p);
        }
        for(auto &e: innerTopEdges){
            outerTopEdges.push_back({e[0]+base, e[1]+base});
        }
        triangulateio topFacet;
        SurfaceMesh topSurfaceMesh;
        generateTRIANGULATEIOWithEdges(outerTopEdgeNodes, outerTopEdges, topHoles, 0, topFacet);
        topSurfaceMesh.projectTRIANGULATEIO(topFacet, PROJECTION_TYPE::XY_PLANE, upperBound[2]);
        deleteTRIANGULATEIOAllocatedArrays(topFacet);

        std::vector<std::array<double,3>> bottomHoles;
        std::array<double,3> bottomHole={0,0, 0};
        double bottomFactor = 1.0/ double(innerBottomEdges.size());
        for(auto &n: innerBottomEdgeNodes){
            bottomHole[0]+=bottomFactor*n[0];
            bottomHole[1]+=bottomFactor*n[1];
        }
        bottomHoles.push_back(bottomHole);

        base=outerBottomEdgeNodes.size();
        for(auto &p: innerTopEdgeNodes){
            outerBottomEdgeNodes.push_back(p);
        }
        for(auto &e: innerTopEdges){
            outerBottomEdges.push_back({e[0]+base, e[1]+base});
        }
        triangulateio bottomFacet;
        SurfaceMesh bottomSurfaceMesh;
        generateTRIANGULATEIOWithEdges(outerBottomEdgeNodes, outerBottomEdges, bottomHoles, 0, bottomFacet);
        bottomSurfaceMesh.projectTRIANGULATEIO(bottomFacet, PROJECTION_TYPE::XY_PLANE, lowerBound[2]);
        deleteTRIANGULATEIOAllocatedArrays(bottomFacet);

        resultingMesh.addSubSurfaceMesh(topSurfaceMesh);
        resultingMesh.addSubSurfaceMesh(bottomSurfaceMesh);
        
    }

    void generateBoundingBoxSurfaceMesh(
        const Vector3D &lowerBound,
        const Vector3D &upperBound,
        double size,
        SurfaceMesh &resultingSurfaceMesh
    ){
        std::vector<std::array<int, 2>> bottomEdges;
        std::vector<std::array<double, 3>> bottomEdgeNodes;

        SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
        std::vector<std::array<double,3>> holes;
        generateRectangleEdges({upperBound[0],upperBound[1]},{lowerBound[0], lowerBound[1]}, size, bottomEdgeNodes, bottomEdges);

        triangulateio triTopBottom;
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triTopBottom);
        faceTop.projectTRIANGULATEIO(triTopBottom, PROJECTION_TYPE::XY_PLANE, upperBound[2]);
        faceBottom.projectTRIANGULATEIO(triTopBottom, PROJECTION_TYPE::XY_PLANE, lowerBound[2]);
        deleteTRIANGULATEIOAllocatedArrays(triTopBottom);

        bottomEdgeNodes.clear();
        bottomEdges.clear();
        holes.clear();
        generateRectangleEdges({upperBound[1], upperBound[2]}, {lowerBound[1],lowerBound[2]}, size, bottomEdgeNodes, bottomEdges);
        triangulateio triFrontBack;	
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triFrontBack);
        faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, upperBound[0]);
        faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, lowerBound[0]);
        deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

        bottomEdgeNodes.clear();
        bottomEdges.clear();
        holes.clear();
        generateRectangleEdges({upperBound[2], upperBound[0]}, {lowerBound[2],lowerBound[0]}, size, bottomEdgeNodes, bottomEdges);
        triangulateio triLeftRight;	
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triLeftRight);
        faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, upperBound[1]);
        faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, lowerBound[1]);
        deleteTRIANGULATEIOAllocatedArrays(triLeftRight);

        double eps = 0.01*size;
        resultingSurfaceMesh.mergeSurfaceMesh(faceRight,eps);
        resultingSurfaceMesh.mergeSurfaceMesh(faceLeft, eps);
        resultingSurfaceMesh.mergeSurfaceMesh(faceFront, eps);
        resultingSurfaceMesh.mergeSurfaceMesh(faceBack, eps);
        resultingSurfaceMesh.mergeSurfaceMesh(faceTop, eps);	
        resultingSurfaceMesh.mergeSurfaceMesh(faceBottom, eps);        

        resultingSurfaceMesh.rebuildIndices();
    }

    void removeTopBottomTrianglesInSurfaceMesh(
        SurfaceMesh &surfaceMesh,
        double top,
        double bottom,
        std::vector<std::array<double, 3>> &topEdgeNodes,
        std::vector<std::array<double, 3>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges,
        std::vector<std::array<int, 2>> &bottomEdges,
        double eps
    ){

        surfaceMesh.rebuildTriangleAdjacency();
        std::unordered_map<int, int> topNodeMap;
        std::unordered_map<int, int> bottomNodeMap;
        int topIndex = 0;
        int bottomIndex = 0;
        auto getTopNode =
        [&surfaceMesh, &topEdgeNodes, &topNodeMap, &topIndex]
        (int index){
            int rst = 0;
            if(topNodeMap.count(index)){
                rst = topNodeMap[index];
            }
            else{
                rst = topIndex++;
                topNodeMap[index]  = rst;
                topEdgeNodes.push_back({surfaceMesh.nodes[index]->pos[0], surfaceMesh.nodes[index]->pos[1], surfaceMesh.nodes[index]->pos[2]});
            }
            return rst;
        };
        auto getBottomNode =
        [&surfaceMesh, &bottomEdgeNodes, &bottomNodeMap, &bottomIndex]
        (int index){
            int rst = 0;
            if(bottomNodeMap.count(index)){
                rst = bottomNodeMap[index];
            }
            else{
                rst = bottomIndex++;
                bottomNodeMap[index]  = rst;
                bottomEdgeNodes.push_back({surfaceMesh.nodes[index]->pos[0], surfaceMesh.nodes[index]->pos[1], surfaceMesh.nodes[index]->pos[2]});
            }
            return rst;
        };
        
        for (auto &n: surfaceMesh.nodes){
            if ((n->pos[2]+eps)>=top){
                n->edit = 1;
            }
            else if((n->pos[2]-eps)<=bottom){
                n->edit = 2;
            }
            else{
                n->edit = 0;
            }
        }
        std::vector<Triangle *> delTriangles;
        for (auto &t: surfaceMesh.triangles){
            int count1 = 0;
            int count2 = 0;
            t->edit = 0;
            for(auto n: t->nodes){
                if (n->edit==1){
                    count1++;
                }
                else if(n->edit==2){
                    count2++;
                }
            }
            if(count1==3){
                t->edit = 1;
                delTriangles.push_back(t);
            }
            else if (count2==3){
                t->edit = 2;
                delTriangles.push_back(t);
            }
        }
        for(auto &t: surfaceMesh.triangles){

            if (t->edit){
                for(int i=0; i<3; i++){
                    Triangle *tt = t->adjacentTriangles[i];
                    if (tt->edit==0){

                        if (t->edit==1){
                            topEdges.push_back({getTopNode(t->nodes[(i+1)%3]->index), getTopNode(t->nodes[(i+2)%3]->index)});
                        }
                        else{
                            bottomEdges.push_back({getBottomNode(t->nodes[(i+1)%3]->index), getBottomNode(t->nodes[(i+2)%3]->index)});						
                        }
                    }
                }			
            }
        }
        surfaceMesh.deleteTriangles(delTriangles);
    }

    void generateBoundingBoxSurfaceMeshWithXYHoles(
        const Vector3D &lowerBound,
        const Vector3D &upperBound,
        double size,
        std::vector<std::array<double, 3>> &topEdgeNodes,
        std::vector<std::array<double, 3>> &bottomEdgeNodes,
        std::vector<std::array<int, 2>> &topEdges, 
        std::vector<std::array<int, 2>> &bottomEdges,
        SurfaceMesh &resultingSurfaceMesh,
        bool hasHoles

    ){
        hasHoles= true;
        SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
        std::array<double, 2> oxymax({upperBound[0], upperBound[1]});
        std::array<double, 2> oxymin({lowerBound[0], lowerBound[1]});

        std::vector<std::array<double,3>> topHoles;
        std::vector<std::array<double,3>> bottomHoles;
        std::vector<std::vector<std::array<double,3>>> topInnerLoopNodes;
        std::vector<std::vector<std::array<int, 2>>> topInnerLoopEdges;
        std::vector<std::vector<std::array<double,3>>> bottomInnerLoopNodes;
        std::vector<std::vector<std::array<int, 2>>> bottomInnerLoopEdges;
        findHoles(topEdgeNodes, topEdges, topInnerLoopNodes, topInnerLoopEdges, topHoles);
        findHoles(bottomEdgeNodes, bottomEdges, bottomInnerLoopNodes, bottomInnerLoopEdges, bottomHoles);
        int nTop= topEdgeNodes.size();
        int nBottom= bottomEdgeNodes.size();

        // std::array<double,3> topHole={0,0,0};
        // double topFactor = 1.0/ double(topEdgeNodes.size());
        // for(auto &n: topEdgeNodes){
        //     topHole[0]+=topFactor*n[0];
        //     topHole[1]+=topFactor*n[1];
        // }
        // topHoles.push_back(topHole);


        // std::array<double,3> bottomHole={0,0,0};
        // double bottomFactor = 1.0/ double(bottomEdges.size());
        // for(auto &n: bottomEdgeNodes){
        //     bottomHole[0]+=bottomFactor*n[0];
        //     bottomHole[1]+=bottomFactor*n[1];
        // }
        // bottomHoles.push_back(bottomHole);

        generateRectangleEdges({upperBound[0],upperBound[1]},{lowerBound[0], lowerBound[1]}, size, bottomEdgeNodes, bottomEdges);
        generateRectangleEdges({upperBound[0],upperBound[1]},{lowerBound[0], lowerBound[1]}, size, topEdgeNodes, topEdges);
        
        triangulateio triBottom;
        triangulateio triTop;	
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, bottomHoles, size*size/2, triBottom);
        generateTRIANGULATEIOWithEdges(topEdgeNodes, topEdges, topHoles, size*size/2, triTop);

        faceBottom.projectTRIANGULATEIO(triBottom, PROJECTION_TYPE::XY_PLANE, lowerBound[2]);
        faceTop.projectTRIANGULATEIO(triTop, PROJECTION_TYPE::XY_PLANE, upperBound[2]);
        deleteTRIANGULATEIOAllocatedArrays(triBottom);
        deleteTRIANGULATEIOAllocatedArrays(triTop);
        for(int i=0; i<nTop; i++){
            faceTop.nodes[i]->pos[2] = topEdgeNodes[i][2];
        }
        faceTop.smooth(3);
        for(int i=0; i<nBottom; i++){
            faceBottom.nodes[i]->pos[2] = bottomEdgeNodes[i][2];
        }
        faceBottom.smooth(3);
        std::vector<std::array<double,3>> holes;
        bottomEdgeNodes.clear();
        bottomEdges.clear();
        generateRectangleEdges({upperBound[1], upperBound[2]}, {lowerBound[1],lowerBound[2]}, size, bottomEdgeNodes, bottomEdges);
        triangulateio triFrontBack;	
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triFrontBack);
        faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, upperBound[0]);
        faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, lowerBound[0]);
        deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

        bottomEdgeNodes.clear();
        bottomEdges.clear();
        generateRectangleEdges({upperBound[2], upperBound[0]}, {lowerBound[2],lowerBound[0]}, size, bottomEdgeNodes, bottomEdges);
        triangulateio triLeftRight;	
        generateTRIANGULATEIOWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triLeftRight);
        faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, upperBound[1]);
        faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, lowerBound[1]);
        deleteTRIANGULATEIOAllocatedArrays(triLeftRight);


        auto addNodeLabel=[](auto &mesh, int label){
            for(auto n: mesh.nodes){
                n->label = label;
            }
        };
  
        addNodeLabel(faceTop, 1);
        addNodeLabel(faceBottom, 1);
        addNodeLabel(faceFront, 1);
        addNodeLabel(faceBack, 1);
        addNodeLabel(faceRight, 1);
        addNodeLabel(faceLeft, 1);
        resultingSurfaceMesh.mergeSurfaceMesh(faceRight, 1e-8);
        resultingSurfaceMesh.mergeSurfaceMesh(faceLeft, 1e-8);
        resultingSurfaceMesh.mergeSurfaceMesh(faceFront, 1e-8);
        resultingSurfaceMesh.mergeSurfaceMesh(faceBack, 1e-8);
        resultingSurfaceMesh.mergeSurfaceMesh(faceTop, 1e-8);
        resultingSurfaceMesh.mergeSurfaceMesh(faceBottom, 1e-8);

        if (topInnerLoopNodes.size()){
            hasHoles = false;
            std::vector<std::array<double,3>> tmpHoles;
            for(int i = 0; i < topInnerLoopNodes.size(); i++){
                triangulateio tmp;
                generateTRIANGULATEIOWithEdges(topInnerLoopNodes[i], topInnerLoopEdges[i], tmpHoles, size*size/2, tmp);
                SurfaceMesh tmpSurfaceMesh;
                tmpSurfaceMesh.projectTRIANGULATEIO(tmp, PROJECTION_TYPE::XY_PLANE, upperBound[2]);
                
                for(int j=0; j<topInnerLoopNodes[i].size(); j++){
                    tmpSurfaceMesh.nodes[j]->pos[2] = topInnerLoopNodes[i][j][2];
                }
                tmpSurfaceMesh.smooth(3);
                resultingSurfaceMesh.mergeSurfaceMesh(tmpSurfaceMesh, 1e-8);
            }

        }
        if (bottomInnerLoopNodes.size()){
            hasHoles = false;
            std::vector<std::array<double,3>> tmpHoles;
            for(int i = 0; i < bottomInnerLoopNodes.size(); i++){
                triangulateio tmp;
                generateTRIANGULATEIOWithEdges(bottomInnerLoopNodes[i], bottomInnerLoopEdges[i], tmpHoles, size*size/2, tmp);
                SurfaceMesh tmpSurfaceMesh;
                tmpSurfaceMesh.projectTRIANGULATEIO(tmp, PROJECTION_TYPE::XY_PLANE, lowerBound[2]);
                
                for(int j=0; j<bottomInnerLoopNodes[i].size(); j++){
                    tmpSurfaceMesh.nodes[j]->pos[2] = bottomInnerLoopNodes[i][j][2];
                }
                tmpSurfaceMesh.smooth(3);
                resultingSurfaceMesh.mergeSurfaceMesh(tmpSurfaceMesh, 1e-8);
            }
        }
    }

    void removeIntersectionElements(
        Mesh &originalMesh,
        Mesh &expandedAtomisticMesh,
        SurfaceMesh &interfaceSurfaceMesh,
        int outerLayers
    ){
        expandedAtomisticMesh.readyForSpatialSearch();

        for(auto n: originalMesh.nodes){
            n->edit = 0;
        }

        for(auto e: originalMesh.tetrahedrons){
            e->edit = 0;
            if (e->label==ATOMISTIC_TET){
                e->edit = 1;
                continue;
            }
            for(auto n: e->nodes){
                if (n->edit==1){
                    e->edit = 1;
                    break;
                }


                if (expandedAtomisticMesh.checkTetrahedronContain(n->pos)){
                    e->edit = 1;
                    n->edit = 1;
                }

                if(e->edit==1) break;
            }

            if(e->edit == 0){
                if(expandedAtomisticMesh.checkTetrahedronIntersection(e)){
                    e->edit = 1;
                }
            }
        }

        std::vector<Tetrahedron *> delTets;
        for(auto e: originalMesh.tetrahedrons){
            if(e->edit==1){
                delTets.push_back(e);
            }
        }

        if (outerLayers>0){
            std::vector<Tetrahedron *> removeTets;
            for(auto n: originalMesh.nodes){
                n->edit = 0;
            }

            for(auto e: originalMesh.tetrahedrons){
                if(e->edit==1){
                    removeTets.push_back(e);
                }
            }
            for(auto e: removeTets){
                for(auto n: e->nodes){
                    n->edit = 1;
                }
            }


            for(int i=0; i<outerLayers; i++){
                std::vector<Tetrahedron *> freshTets;
                for(auto e: originalMesh.tetrahedrons){
                    if (e->edit == 0){
                        for(auto n: e->nodes){
                            if (n->edit == 1){
                                e->edit = 1;
                                freshTets.push_back(e);
                                break;
                            }
                        }
                    }
                }

                for(auto e: freshTets){
                    for(auto n: e->nodes){
                        n->edit = 1;
                    }
                }

                if(freshTets.empty()){
                    break;
                }
            }
        }

        extractBorderingSurfaceMesh(delTets, interfaceSurfaceMesh);
        originalMesh.deleteTetrahedrons(delTets);

    }


    void extractBorderingSurfaceMesh(
        std::vector<Tetrahedron *>&tets, 
        SurfaceMesh &aSurface
    ){
        std::set<Node *> nodeSet; 
        std::unordered_set<SubTriangle, SubTriangleHasher, SubTriangleEqual> subTriangleSet;
        for(auto e: tets){
            for(int i=0; i<4; i++){
                SubTriangle keyFacet = e->getSubTriangle(i);
                if (subTriangleSet.count(keyFacet)){
                    subTriangleSet.erase(keyFacet);
                }
                else{
                    subTriangleSet.insert(keyFacet);
                }
            }
        }
        std::unordered_map<Node*, Node*> oldNewNodes;
        auto getNode
        =
        [&oldNewNodes]
        (Node* key){
            Node *rst;
            if(oldNewNodes.find(key)!=oldNewNodes.end()){
                rst = oldNewNodes[key];
            }
            else{
                rst = new Node(key->pos);
                rst->label= key->label;
                oldNewNodes[key] = rst;
            }
            return rst;
        };

        for(auto f: subTriangleSet){
            aSurface.addTriangle(getNode(f.forms[0]), getNode(f.forms[1]), getNode(f.forms[2]));
            
        }
        for(auto kv: oldNewNodes){
            aSurface.nodes.push_back(kv.second);
        }

        aSurface.rebuildIndices();
    }
}