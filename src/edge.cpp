#include "edge.h"
#include "node.h"

void Edge::sortNodeIndices(){
        int i0 = sNodes[0]->index;
        int i1 = sNodes[1]->index;
        if (i0<i1){
            orderedNodeIndices[0] = i0;
            orderedNodeIndices[1] = i1;
        }
        else{
            orderedNodeIndices[0] = i1;
            orderedNodeIndices[1] = i0;
        }    
}