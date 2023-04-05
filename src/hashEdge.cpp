#include "hashEdge.h"



void HashEdgeTable::insert(const Edge &aEdge){
    columns[key(aEdge)].push_back(aEdge);
}

int HashEdgeTable::key(const Edge &aEdge){
    return 3*aEdge.orderedNodeIndices[0]+7*aEdge.orderedNodeIndices[1];
}

bool HashEdgeTable::search(Edge &keyEdge, std::vector<Edge> &goalEdge){
    int theKey = key(keyEdge);
    bool rst = false;
    for(auto &e: columns[theKey]){
        if(e.orderedNodeIndices[0] == keyEdge.orderedNodeIndices[0]
        && e.orderedNodeIndices[1] == keyEdge.orderedNodeIndices[1]){
            goalEdge.push_back(e);
            rst=true;
        }
    }
    return rst;
}

bool HashEdgeTable::exist(const Edge &keyEdge){
    int theKey = key(keyEdge);
    bool rst = false;
    for(auto &e: columns[theKey]){
        if(e.orderedNodeIndices[0] == keyEdge.orderedNodeIndices[0]
        && e.orderedNodeIndices[1] == keyEdge.orderedNodeIndices[1]){
            rst=true;
            break;
        }
    }
    return rst;
}




bool HashEdgeTable::remove(int key, const Edge& keyEdge){
    bool rst = false;
    for(int i=0; i<columns[key].size(); i++){
        Edge& e= columns[key][i];
        if(e.orderedNodeIndices[0] == keyEdge.orderedNodeIndices[0]
        && e.orderedNodeIndices[1] == keyEdge.orderedNodeIndices[1]){
            columns[key].erase(columns[key].begin()+i);
            i--;
            rst = true;
        }
    }
    return rst;
}

bool HashEdgeTable::remove(const Edge &keyEdge){
    bool rst = false;
    int theKey = key(keyEdge);
    for(int i=0; i<columns[theKey].size(); i++){
        Edge& e= columns[theKey][i];
        if(e.orderedNodeIndices[0] == keyEdge.orderedNodeIndices[0]
        && e.orderedNodeIndices[1] == keyEdge.orderedNodeIndices[1]){
            columns[theKey].erase(columns[theKey].begin()+i);
            i--;
            rst = true;
        }
    }
    return rst;
}