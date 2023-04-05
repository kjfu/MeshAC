#include "surfaceMesh.h"
#include "common.h"
#include <fstream>
#include "hashEdge.h"
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "kdtree.h"

void SurfaceMesh::exportVTK(const std::string &filePath){
    std::ofstream file(filePath);
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "vtk output"                 << std::endl;
    file << "ASCII"                      << std::endl;
    file << "DATASET UNSTRUCTURED_GRID"  << std::endl;
    file << "POINTS " << nodes.size() << " double" << std::endl;
    for(auto n: nodes){
        file <<n->pos[0] <<  "  " << n->pos[1] << "  " << n->pos[2] << std::endl;
    }

    file << "CELLS " << triangles.size() << "  " << 4* triangles.size() << std::endl;
    for(auto tri: triangles){
        file << "3 " << tri->nodes[0]->index
        << "  " << tri->nodes[1]->index
        << "  " << tri->nodes[2]->index
        <<std::endl;
    }
    file << "CELL_TYPES " << triangles.size() << std::endl;
    for(auto tri: triangles){
        file << VTK_ELEMENT_TYPE::VTK_TRIANGLE << std::endl;
    }

    file.close();
}






void SurfaceMesh::rebuildIndices(){
    for(int i=0; i<nodes.size(); i++){
        nodes[i]->index = i;
    }

    for(int i=0; i<triangles.size(); i++){
        triangles[i]->index = i;
    }    
}


void SurfaceMesh::exportMESH(const std::string &filePath){
    std::ofstream file(filePath);
    file << "MeshVersionFormatted 2\n" << std::endl;
    file << "Dimension 3\n" << std::endl;
    file << "Vertices\n" << nodes.size() << std::endl;
    for(auto& n: nodes){
        file << n->pos[0] << " " << n->pos[1] << " " << n->pos[2] << " " << n->label << std::endl;
    }
    
    file << "Triangles\n" << triangles.size() << std::endl;
    for(auto& t: triangles){
        file << t->nodes[0]->index+1
            << "  " << t->nodes[1]->index+1
            << "  " << t->nodes[2]->index+1
            << "  " << t->label << std::endl;
    }

    file.close();
}

void SurfaceMesh::loadMESH(const std::string &filePath){
    std::ifstream file(filePath);
    if(file.is_open()){
        while(file){
            std::string line;
            std::getline(file, line);
            std::stringstream lineStream(line);
            std::string keyString;
            lineStream >> keyString;
            if (keyString=="Vertices"){
                std::string numLine;
                std::getline(file, numLine);
                int nv = std::stoi(numLine);
                for(int i=0; i<nv; i++){
                    std::string vertexLine;
                    std::getline(file, vertexLine);
                    std::stringstream vertexStream(vertexLine);
                    double x, y, z;
                    int label;
                    vertexStream >> x >> y >> z >> label;
                    Node *n = new Node(x, y, z);
                    n->label = label;
                    nodes.push_back(n);
                }
            }
            else if (keyString == "Triangles"){
                std::string numLine;
                std::getline(file, numLine);
                int nt = std::stoi(numLine);
                for(int i=0; i<nt; i++){
                    std::string triangleLine;
                    std::getline(file, triangleLine);
                    std::stringstream triangleStream(triangleLine);
                    int n0, n1, n2;
                    int label;
                    triangleStream >> n0 >> n1 >> n2 >> label;
                    TriangleElement *e = new TriangleElement(nodes[n0-1], nodes[n1-1], nodes[n2-1]);
                    e->label = label;
                    triangles.push_back(e);
                }
            }
        }
        rebuildIndices();
        file.close();
    }
}



void SurfaceMesh::rebuildTriangleAdjacency(){
    HashEdgeTable anEdgeTable;
    
    for (auto &t: triangles){
        for(int i=0; i<3; i++){
            t->adjacentTriangles[i].clear();
            Edge e = t->edge(i);
            anEdgeTable.insert(e);
        }
    }

    for(auto &t: triangles){
        for(int i=0; i<3; i++){
            if (t->adjacentTriangles[i].empty()){
                std::vector<Edge> edges;
                Edge keyEdge = t->edge(i);
                if(anEdgeTable.search(keyEdge, edges)){
                    for(int j=0; j<edges.size(); j++){
                        TriangleElement* currentTriangle = edges[j].tri;
                        for(int k=0; k<edges.size(); k++){
                            Edge goalEdge = edges[k];
                            if(goalEdge.tri != currentTriangle){
                                goalEdge.tri->adjacentTriangles[goalEdge.localIndexOfTriangle].push_back(currentTriangle);
                            }
                        }
                    }
                }
                
            } 
        }
    }

}


// void SurfaceMesh::deleteTriangles(std::vector<TriangleElement*> delTriangles){
//     for(auto t: triangles){
//         t->edit = 0;
//         for(auto n:t->nodes){
//             n->edit = 0;
//         }
//     }

//     for(auto t: delTriangles){
//         t->edit=-1;
//         for(auto n: t->nodes){
//             n->edit = -1;
//         }
//     }
//     for(auto t:triangles){
//         if(t->edit==0){
//             for(int i=0; i<3; i++){
//                 for(int j=0; j<t->adjacentTriangles[i].size(); j++){
//                     if(t->adjacentTriangles[i][j]->edit==-1){
//                         t->adjacentTriangles[i].erase(t->adjacentTriangles[i].begin()+j);
//                         j--;
//                     }
//                 }
//             }
//             for(auto n:t->nodes){
//                 n->edit = 0;
//             }
//         }
//     }

//     for(int i=0; i<nodes.size(); i++){
//         if(nodes[i]->edit==-1){
//             delete nodes[i];
//             nodes.erase(nodes.begin()+i);
//             i--;
//         }
//     }

//     for(int i=0; i<triangles.size(); i++){
//         if(triangles[i]->edit==-1){
//             delete triangles[i];
//             triangles.erase(triangles.begin()+i);
//             i--;
//         }
//     }

//     rebuildIndices();
// }




void SurfaceMesh::loadPLY(const std::string &filePath){
    std::ifstream file(filePath);
    if (file.is_open()){
        while(file){
            std::string line;
            std::getline(file, line);
            std::stringstream lineStream(line);
            std::string key0, key1;
            lineStream >> key0 >> key1;
            int nv, nf;
            if (key0=="element" && key1=="vertex"){
                lineStream >> nv;
            }
            else if(key0=="element" && key1=="face"){
                lineStream >> nf;
            }
            else if (key0=="end_header"){
                for(int i=0; i<nv; i++){
                    std::string vertexLine;
                    std::getline(file, vertexLine);
                    std::stringstream vertexStream(vertexLine);
                    double x, y, z;
                    vertexStream >> x >> y >> z;
                    Node* n = new Node(x, y, z);
                    nodes.push_back(n);
                }
                
                for(int i=0; i<nf; i++){
                    std::string faceLine;
                    std::getline(file, faceLine);
                    std::stringstream faceStream(faceLine);
                    int tmp, n0, n1, n2;
                    faceStream >> tmp >> n0 >> n1 >> n2;
                    TriangleElement *t = new TriangleElement(nodes[n0], nodes[n1], nodes[n2]);
                    triangles.push_back(t);
                }
            }
        }
        rebuildIndices();
        file.close();
    }
}

void SurfaceMesh::exportPLY(const std::string &filePath){
    std::ofstream file(filePath);
    file << "ply" << std::endl;
    file << "format ascii 1.0" << std::endl;
    file << "comment VCGLIB generated" << std::endl;
    file << "element vertex " << nodes.size() << std::endl;
    file << "property float x" << std::endl;
    file << "property float y" << std::endl;
    file << "property float z" << std::endl;
    file << "element face "<< triangles.size() << std::endl;
    file << "property list uchar int vertex_indices" << std::endl;
    file << "end_header" << std::endl;
    for(int i=0; i<nodes.size(); i++){
        file << nodes[i]->pos[0]
        << " " << nodes[i]->pos[1] 
        << " " << nodes[i]->pos[2] << std::endl;
    }

    for(int i=0; i<triangles.size(); i++){
        file << 3 
        << " " << triangles[i]->nodes[0]->index
        << " " << triangles[i]->nodes[1]->index
        << " " << triangles[i]->nodes[2]->index << std::endl;
    }    
}




void SurfaceMesh::loadVTK(const std::string &filePath){
    std::ifstream file(filePath);
    if (file.is_open()){
        while(file){
            std::string line;
            std::getline(file, line);
            std::stringstream lineStream(line);
            std::string keyString;
            lineStream >> keyString;
            if (keyString=="POINTS"){
                int nv;
                lineStream >> nv;
                int get=0;
                while (get<nv)
                {
                    std::string subLine;
                    std::getline(file, subLine);
                    std::stringstream subStream(subLine);
                    while(subStream && get<nv) {
                        
                        std::string  x, y, z;
                        
                        subStream >> x >> y >> z;
                        if(z.size()>0){
                            Node *node = new Node(std::stod(x), std::stod(y),  std::stod(z));
                            nodes.push_back(node);
                            get++;
                        } 
                    }
                }
            }
            else if (keyString=="CELLS"){
                int ne;
                lineStream >> ne;
                int get=0;
                while (get<ne){
                    std::string subLine;
                    std::getline(file, subLine);
                    std::stringstream subStream(subLine);
                    while(subStream && get<ne) {
                        std::string tmp, n0, n1, n2;
                        subStream >> tmp >> n0 >> n1 >> n2 ;
                        if(n2.size()>0){
                            TriangleElement *t = new TriangleElement(nodes[std::stoi(n0)], nodes[std::stoi(n1)], nodes[std::stoi(n2)]);
                            triangles.push_back(t);
                            get++; 
                        } 
                    }
                }
                
            }


        }
        rebuildIndices();
        file.close();
    }    
}

void SurfaceMesh::addSubSurfaceMesh(SurfaceMesh &another){
    int base = this->nodes.size();
    for(auto n: another.nodes){
        Node *newNode = new Node(n->pos);
        newNode->label= n->label;
        
        this->nodes.push_back(newNode);
    }

    for(auto e: another.triangles){
        std::vector<Node *> nodesToAdd;
        for(auto n: e->nodes){
            nodesToAdd.push_back(this->nodes[base+n->index]);
        
        }

        TriangleElement *tri = new TriangleElement(nodesToAdd[0], nodesToAdd[1], nodesToAdd[2]);
        tri->label = e->label;
        this->triangles.push_back(tri);
    }

    this->rebuildIndices();

}

void SurfaceMesh::exportSU2(const std::string &filePath){
    std::ofstream file(filePath);
    
    file << "NDIME= 3"  << std::endl;

    file << "NELEM= 0" << std::endl; 
    file << "NPOIN= " << nodes.size() << std::endl;
    for(int i=0; i<nodes.size(); i++){
        file << nodes[i]->pos[0] 
        << "   " << nodes[i]->pos[1] 
        << "   " << nodes[i]->pos[2] 
        << "   " << i << std::endl;
    }


    std::unordered_map<int, std::vector<TriangleElement *>> zones;
    for(auto e: triangles){
        zones[e->label].push_back(e);
    }
    file << "NMARK= "<< zoneNames.size() << std::endl;
    for(auto kv:zoneNames){
        file << "MARKER_TAG= " << zoneNames[kv.first]<< std::endl;
        file << "MARKER_ELEMS= " << zones[kv.first].size() << std::endl;
        for(auto e: zones[kv.first]){
            file << VTK_ELEMENT_TYPE::VTK_TRIANGLE 
            << "  " << e->nodes[0]->index
            << "  " << e->nodes[1]->index
            << "  " << e->nodes[2]->index
            << std::endl;
        }

    }


}

void SurfaceMesh::exportTETGENIO(tetgenio &out){
    out.numberofpoints = this->nodes.size();
    out.firstnumber = 0;

    out.pointlist = new REAL[out.numberofpoints*3];
    out.pointmarkerlist = new int[out.numberofpoints];
    for(int i=0; i<out.numberofpoints; i++){
        out.pointlist[3*i] = nodes[i]->pos[0];
        out.pointlist[3*i+1] = nodes[i]->pos[1];
        out.pointlist[3*i+2] = nodes[i]->pos[2];
        out.pointmarkerlist[i] = nodes[i]->label;
    }

    out.numberoffacets = triangles.size();
    out.facetlist = new tetgenio::facet[out.numberoffacets];
    out.facetmarkerlist = new int[out.numberoffacets];
    for(int i=0; i<out.numberoffacets; i++){
        tetgenio::facet *f = &out.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = nullptr;
        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = triangles[i]->nodes[0]->index;
        p->vertexlist[1] = triangles[i]->nodes[1]->index;
        p->vertexlist[2] = triangles[i]->nodes[2]->index;

        out.facetmarkerlist[i] = triangles[i]->label;
    }
}


void  SurfaceMesh::checkNonManifoldEdges(std::vector<Edge> &nonManifoldEdges){
    HashEdgeTable anEdgeTable;
    
    for (auto &t: triangles){
        for(int i=0; i<3; i++){
            t->adjacentTriangles[i].clear();
            Edge e = t->edge(i);
            anEdgeTable.insert(e);
        }
    }


    for(auto &t: triangles){
        for(int i=0; i<3; i++){
            if (t->adjacentTriangles[i].empty()){
                std::vector<Edge> edges;
                Edge keyEdge = t->edge(i);
                if(anEdgeTable.search(keyEdge, edges)){
                    for(int j=0; j<edges.size(); j++){
                        TriangleElement* currentTriangle = edges[j].tri;
                        for(int k=0; k<edges.size(); k++){
                            Edge goalEdge = edges[k];
                            if(goalEdge.tri != currentTriangle){
                                goalEdge.tri->adjacentTriangles[goalEdge.localIndexOfTriangle].push_back(currentTriangle);
                            }
                        }
                    }
                    if(edges.size()>2){
                        nonManifoldEdges.push_back(keyEdge);
                    }
                }
            } 
        }
    }

}



void SurfaceMesh::deleteTriangles(std::vector<TriangleElement *> &toDelTriangles){
    for(auto n: nodes){
        n->edit = 1;
    }

    for(auto e:triangles){
        e->edit = 0;
    }

    for(auto e: toDelTriangles){
        e->edit=1;
    }

    for(auto e:triangles){
        if(e->edit==0){
            for(auto n:e->nodes){
                n->edit=0;
            }
        }
    }

    for(int i=0; i<triangles.size(); i++){
        if(triangles[i]->edit){
            delete triangles[i];
            triangles.erase(triangles.begin()+i);
            i--;
        }
    }

    for(int i=0; i<nodes.size(); i++){
        if(nodes[i]->edit){
            delete nodes[i];
            nodes.erase(nodes.begin()+i);
            i--;
        }
    }

    rebuildIndices();
}


void SurfaceMesh::mergeSurfaceMesh(SurfaceMesh &another, double tolerance){
    struct kdtree *kd = kd_create(3);
    
    for(auto n: nodes){
        kd_insert(kd, n->pos.data(), n);
    }

    std::unordered_map<Node*, Node*> nodeMap;
    std::vector<Node *> addNodes;
    int numMergeNodes=0;
    for(auto n: another.nodes){
        kdres *search_rst =  kd_nearest_range(kd, n->pos.data(),  tolerance);
        if (kd_res_size(search_rst)){
            numMergeNodes++;
            void *data = kd_res_item_data(search_rst);
            Node *nn = static_cast<Node*>(data);
            nodeMap[n] = nn;

        }
        else{
            Node *nn = new Node(*n);
            nn->label = n->label;
            addNodes.push_back(nn);
            nodeMap[n]= nn;
        }
        kd_res_free(search_rst);
    }
    kd_free(kd);

    nodes.insert(nodes.end(), addNodes.begin(), addNodes.end());
    for(auto e: another.triangles){
        TriangleElement *t = new TriangleElement(
        nodeMap[e->nodes[0]]
        , nodeMap[e->nodes[1]]
        , nodeMap[e->nodes[2]]);
        t->label = e->label;
        triangles.push_back(t);
    }

    rebuildIndices();
}



void SurfaceMesh::projectTRIANGULATEIO(triangulateio &in, PROJECTION_TYPE projectionType, double offset){
    int u, v, w;
    if (projectionType == PROJECTION_TYPE::XY_PLANE){
        u=0; v=1; w=2;
    }
    else if (projectionType == PROJECTION_TYPE::YZ_PLANE){
        u=1; v=2; w=0;
    }
    else if(projectionType == PROJECTION_TYPE::ZX_PLANE){
        u=2; v=0; w=1;
    }

    for(int i=0; i<in.numberofpoints; i++){
        std::array<double, 3> pos;
        pos[u] = in.pointlist[2*i];
        pos[v] = in.pointlist[2*i+1];
        pos[w] = offset;
        Node *n = new Node(pos.data());
        nodes.push_back(n);
    }

    for(int i=0; i<in.numberoftriangles; i++){
        TriangleElement *t = new TriangleElement(nodes[in.trianglelist[3*i]-1], nodes[in.trianglelist[3*i+1]-1], nodes[in.trianglelist[3*i+2]-1]);
        triangles.push_back(t);
    }

    rebuildIndices();

}


void SurfaceMesh::cloneSurfaceMesh(SurfaceMesh &another){
    std::unordered_map<Node*, Node*> nodeMap;
    for(auto n: another.nodes){
        Node *nn = new Node(n->pos);
        nn->label = n->label;
        nodeMap[n] == nn;
        nodes.push_back(nn);
    }

    for(auto e: another.triangles){
        TriangleElement *ee = new TriangleElement(nodeMap[e->nodes[0]], nodeMap[e->nodes[1]], nodeMap[e->nodes[2]]);
        ee->label = e->label;
        triangles.push_back(ee);
    }

    rebuildIndices();
}

void SurfaceMesh::loadTETGENIO(tetgenio &in){

    std::unordered_map<int, Node*> indexNodeMap;
    auto getNode
    =
    [&in, &indexNodeMap]
    (int index){
        if (indexNodeMap.find(index) == indexNodeMap.end()){
            int realIndex  = index - in.firstnumber;
            indexNodeMap[index] = new Node(in.pointlist[3*realIndex], in.pointlist[3*realIndex+1],in.pointlist[3*realIndex+2]);
        }

        return indexNodeMap[index];
    };


    for(int i=0; i<in.numberoftrifaces; i++){
        // if (in.trifacemarkerlist[i]==1){
            TriangleElement *e = new TriangleElement(getNode(in.trifacelist[3*i]), getNode(in.trifacelist[3*i+1]), getNode(in.trifacelist[3*i+2]));
            triangles.push_back(e);
        // }
    }

    for(auto &kv:indexNodeMap){
        nodes.push_back(kv.second);
    }
    
    rebuildIndices();
}


void SurfaceMesh::estimateSizing(){
    maxSizing = std::numeric_limits<double>::min();
    minSizing = std::numeric_limits<double>::max();
    for(auto n: nodes){
        n->sizing = std::numeric_limits<double>::max();
    }
    
    for(auto e: triangles){
        for(int i=0; i<3; i++){
            for(int j=i+1; j<3; j++){
                Node *n0 = e->nodes[i];
                Node *n1 = e->nodes[j];
                Vector3D dist = n0->pos - n1->pos;
                double d = dist.norm();
                if(d<n0->sizing){
                    n0->sizing = d;
                }

                if (d<n1->sizing){
                    n1->sizing = d;
                }

                if (d>maxSizing){
                    maxSizing = d;
                }

                if (d<minSizing){
                    minSizing = d;
                }
            }
        }
    }
}