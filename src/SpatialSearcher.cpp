/*
 * @Author: Kejie Fu
 * @Date: 2022-01-21 17:11:50
 * @LastEditTime: 2023-03-28 16:39:31
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/SpatialSearcher.cpp
 */
#include "SpatialSearcher.h"
#include "mesh.h"


void SpatialSearcher::buildAABBTree(){

    for (auto tet: mesh->tetrahedrons){
        insertTetrahedron(tet);
    }
}


unsigned int SpatialSearcher::getId(){
    unsigned int rst = 0;
    if (trashIds.empty()){
        rst = currentId++;
    }
    else{
        rst = trashIds.top();
        trashIds.pop();
    }
    return rst;
}

void SpatialSearcher:: searchTetrahedronContain(Vector3D position, SearchTetrahedronResult &result){
    std::vector<double> pos = position.toSTDVector();
    aTree.insertParticle(std::numeric_limits<unsigned int>::max(), pos, 0);
    std::vector<unsigned int> ids = aTree.query(std::numeric_limits<unsigned int>::max());
    for(auto id: ids){
        std::vector<int> zeroNodeIndices;
        Tetrahedron *tet = ID2TET[id];
        int insideCounts=0;
        for(int i=0; i<4; i++){
            double o3d = orient3d(
                tet->nodes[TetrahedronFacet[i][0]]->pos.data(), 
                tet->nodes[TetrahedronFacet[i][1]]->pos.data(), 
                tet->nodes[TetrahedronFacet[i][2]]->pos.data(), 
                position.data());
            if (o3d>0) break;
            if (o3d==0) zeroNodeIndices.push_back(i);
            insideCounts++;
            result.weights[i] = o3d;
        }

        if (insideCounts==4){
            result.tet = tet;
            if (zeroNodeIndices.size()==0){
                result.positionType = POSITION_TYPE::INSIDE;
            }
            else if (zeroNodeIndices.size()==1){
                result.positionType = POSITION_TYPE::ONFACE;
                result.localIndex = zeroNodeIndices[0];
            }
            else if(zeroNodeIndices.size()==2){
                result.positionType = POSITION_TYPE::ONEDGE;
                result.localIndex = TetrahedronTwoFacetCommonEdge[zeroNodeIndices[0]][zeroNodeIndices[1]];
            }
            else if (zeroNodeIndices.size()==3){
                result.positionType = POSITION_TYPE::ONVERTEX;
                int totalIndex = 6;
                for(auto ii: zeroNodeIndices){
                    totalIndex -= ii;
                }
                result.localIndex = totalIndex;
            }

            break;
        }
    }



    aTree.removeParticle(std::numeric_limits<unsigned int>::max());
}
void SpatialSearcher::updateTetrahedron(Tetrahedron *pTet){
    if (TET2ID.count(pTet)){
        int id = TET2ID[pTet];
        std::vector<double> lowerBound;
        std::vector<double> upperBound;
        pTet->getBoundingBox(lowerBound, upperBound);
        aTree.updateParticle(id, lowerBound, upperBound);
    }

}
void SpatialSearcher::removeTetrahedron(Tetrahedron *pTet){
    if (TET2ID.count(pTet)){
        int id = TET2ID[pTet];
        trashIds.push(id);
        aTree.removeParticle(id);
        TET2ID.erase(pTet);
        ID2TET.erase(id);
    }
}

void SpatialSearcher::insertTetrahedron(Tetrahedron *pTet){
    std::vector<double> lowerBound;
    std::vector<double> upperBound;
    pTet->getBoundingBox(lowerBound, upperBound);
    int id = getId();
    ID2TET[id] = pTet;
    TET2ID[pTet] = id;
    aTree.insertParticle(id, lowerBound, upperBound);
}


