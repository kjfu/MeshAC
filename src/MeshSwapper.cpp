/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:35:28
 * @LastEditTime: 2023-03-28 15:59:28
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/MeshSwapper.cpp
 */
#include "MeshSwapper.h"
#include <unordered_map>
    void MeshSwapper::SwapOptimization(int SubdomainNumber, double minQuality){
        mesh->rebuildTetrahedronsAdjacency();
        for(auto tet: mesh->tetrahedrons){
            tet->calculateQuality();
            tet->edit=0;
        }

        int MaxIter=5;
        int iter=0;
        int nEdgeSwap=0;
        do{
            nEdgeSwap=0;
            iter++;
            for(int i=0; i<mesh->tetrahedrons.size(); i++){
                if (mesh->tetrahedrons[i]->edit==0 
                && mesh->tetrahedrons[i]->label==SubdomainNumber
                && mesh->tetrahedrons[i]->quality<minQuality){
                    VolumeShell aShell;
                    VolumeUmbrella anUmbrella;
                    std::vector<double> qualities;
                    for (int j=0; j<6; j++){
                        if (checkEdgeSwap(mesh->tetrahedrons[i], j, aShell, anUmbrella, qualities)){
                            EdgeSwap(aShell, anUmbrella, qualities);
                            nEdgeSwap+=aShell.tetrahedronEdges.size();
                            break;
                        }
                    }
                }
            }
            
            //remove deleted tetrahedral elements
            if (nEdgeSwap){
                for(int i=0; i<mesh->tetrahedrons.size(); i++){
                    if (mesh->tetrahedrons[i]->edit==1){
                        delete mesh->tetrahedrons[i];
                        mesh->tetrahedrons[i]=mesh->tetrahedrons.back();
                        mesh->tetrahedrons.pop_back();
                        i--;
                    }
                }
            }
        }
        while(iter<MaxIter && nEdgeSwap);
        
    }


    bool MeshSwapper::checkEdgeSwap(Tetrahedron *tet, int iLocal, VolumeShell &aShell, VolumeUmbrella &anUmbrella, std::vector<double> &qualities){
        mesh->getVolumeShell(tet, iLocal, aShell);
        if (!aShell.closed||aShell.tetrahedronEdges.size()==1) return false;
        std::set<int> labels;
        for(auto &sub: aShell.tetrahedronEdges){
            labels.insert(sub.tetrahedron->label);
        }
        if (labels.size()>1) return false;
        double qBefore = 1000;
        std::vector<SubTriangle> subtriangles;
        for(auto e: aShell.tetrahedronEdges){
            subtriangles.push_back(e.tetrahedron->getSubTriangle(TetrahedronEdge[e.iLocal][0]));
            subtriangles.push_back(e.tetrahedron->getSubTriangle(TetrahedronEdge[e.iLocal][1]));
            if (e.tetrahedron->quality<qBefore){
                qBefore = e.tetrahedron->quality;
            }
        }

        for(auto n: aShell.nodes){
            anUmbrella.endNode = n;
            anUmbrella.tetrahedronTriangles.clear();
            qualities.clear();
            bool success = true;
            for(auto tri: subtriangles){
                if (!tri.contained(n)){
                    double q = calculateTetrahedronQualityWith4Points_ISO(tri.forms[0]->pos, tri.forms[1]->pos, tri.forms[2]->pos, n->pos);
                    if (q<qBefore){
                        success = false;
                    }
                    else{
                        anUmbrella.tetrahedronTriangles.push_back(tri);
                        qualities.push_back(q);
                    }
                }
                if(!success) break;
            }
            if(success) return true;
        }

        return false;
    }

    void MeshSwapper::EdgeSwap(VolumeShell &aShell,VolumeUmbrella &anUmbrella, std::vector<double> &qualities){
        std::unordered_map<SubTriangle, SubTriangle, SubTriangleHasher, SubTriangleEqual> outAdjacentTetrahedronTriangles;
        std::unordered_map<SubTriangle, SubTriangle, SubTriangleHasher, SubTriangleEqual> inAdjacentTetrahedronTriangles;
        for(auto &subedge: aShell.tetrahedronEdges){
            for(int i = 0; i < 2; i++){
                Tetrahedron *adj = subedge.tetrahedron->adjacentTetrahedrons[TetrahedronEdge[subedge.iLocal][i]];
                if (adj){
                    SubTriangle subtriangle = subedge.tetrahedron->getSubTriangle(TetrahedronEdge[subedge.iLocal][i]);
                    SubTriangle adjTri = adj->getSubTriangle(adj->getLocalIndex(subtriangle));
                    outAdjacentTetrahedronTriangles[adjTri] = adjTri;
                }
            }
        }
        int label=aShell.tetrahedronEdges[0].tetrahedron->label;
        Node *keyNode = anUmbrella.endNode;
        for(auto &sub: aShell.tetrahedronEdges){
            sub.tetrahedron->edit=1;
        }
        
        for(int i=0; i<anUmbrella.tetrahedronTriangles.size(); i++){
            SubTriangle &subtri = anUmbrella.tetrahedronTriangles[i];
            Tetrahedron *tmp = mesh->addTetrahedron(subtri.forms[0], subtri.forms[1], subtri.forms[2], keyNode);
            for(int j=0; j<4; j++){
                SubTriangle tmpTri = tmp->getSubTriangle(j);
                if (outAdjacentTetrahedronTriangles.count(tmpTri)){
                    SubTriangle &adj=outAdjacentTetrahedronTriangles[tmpTri];
                    adj.tetrahedron->adjacentTetrahedrons[adj.iLocal]  = tmp;
                    tmp->adjacentTetrahedrons[j] = adj.tetrahedron;
                }
                else if (inAdjacentTetrahedronTriangles.count(tmpTri)){
                    SubTriangle &adj=inAdjacentTetrahedronTriangles[tmpTri];
                    adj.tetrahedron->adjacentTetrahedrons[adj.iLocal]  = tmp;
                    tmp->adjacentTetrahedrons[j] = adj.tetrahedron;
                    inAdjacentTetrahedronTriangles.erase(tmpTri);
                }
                else{
                    inAdjacentTetrahedronTriangles[tmpTri] = tmpTri;
                }
            }

            tmp->quality = qualities[i];
            tmp->label = label;
        }   


    }

