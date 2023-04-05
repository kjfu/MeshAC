/*
 * @Author: Kejie Fu
 * @Date: 2023-03-27 19:42:30
 * @LastEditTime: 2023-03-28 15:57:59
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/MeshRefiner.cpp
 */
#include "MeshRefiner.h"
#include "SubEntities.h"
#include "MeshSwapper.h"
#include <iostream>
void MeshRefiner::refine(int firstIndex, std::vector<int> elementIndices, int subdomain, int nodeLabel){
    double minQ=1000;
    double maxQ=-1;
    for(auto tet: mesh->tetrahedrons){
        tet->edit = 0;
        tet->calculateQuality();
        if (tet->label==subdomain){
            if(tet->quality<minQ){
                minQ = tet->quality;
            }
            if(tet->quality>maxQ){
                maxQ = tet->quality;
            }
        }
    }
    std::cout << "Quality before refine:" << minQ << std::endl;

    for(int i=0; i<elementIndices.size(); i++){
        Tetrahedron *tet= mesh->tetrahedrons[elementIndices[i]-firstIndex];
        
        if (tet->edit==0 && tet->label==subdomain){
            Node* n= mesh->addNode(tet->center());
            n->label = nodeLabel;
            for(int j=0; j<4; j++){
                Tetrahedron *tmp= mesh->addTetrahedron(tet->nodes[TetrahedronFacet[j][0]],
                tet->nodes[TetrahedronFacet[j][1]],
                tet->nodes[TetrahedronFacet[j][2]], 
                n);
                tmp->label = tet->label;
            }
            tet->edit = 1;
        }
    }

    for(int i=0; i<mesh->tetrahedrons.size(); i++){
        if(mesh->tetrahedrons[i]->edit==1){
            delete mesh->tetrahedrons[i];
            mesh->tetrahedrons[i]=mesh->tetrahedrons.back();
            mesh->tetrahedrons.pop_back();
            i--;
        }
    }


    MeshSwapper aSwapper(mesh);
    aSwapper.SwapOptimization(subdomain, minQ+0.3*(maxQ-minQ));
    double minQ2 = 1000;
    for(auto tet: mesh->tetrahedrons){
        if (tet->label==1 && tet->quality<minQ2){
            minQ2 = tet->quality;
        }
    }
    std::cout << "Quality after refine:" << minQ2 << std::endl;

}