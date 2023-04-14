#include "TetgenTool.h"

namespace MeshAC{

    void transportNodesToTETGENIO(const std::vector<Node *> &sNodes, tetgenio &out){
        out.numberofpoints = sNodes.size();
        out.pointlist = new double[out.numberofpoints*3];
        out.pointmarkerlist = new int[out.numberofpoints];
        for(int i=0; i<out.numberofpoints; i++){
            for(int j=0; j<3; j++){
                out.pointlist[3*i+j] = sNodes[i]->pos[j]; 
            }
            out.pointmarkerlist[i] = sNodes[i]->label;
        }
    }

    void transportPointsToTETGENIO(const std::vector<Vector3D> &points, int label, tetgenio &out){
        out.numberofpoints = points.size();
        out.pointlist = new double[out.numberofpoints*3];
        out.pointmarkerlist = new int[out.numberofpoints];
        for(int i=0; i<out.numberofpoints; i++){
            for(int j=0; j<3; j++){
                out.pointlist[3*i+j] = points[i][j]; 
            }
            out.pointmarkerlist[i] = label;
        }
    }

    void transportSurfaceMeshToTETGENIO(
        SurfaceMesh &surfaceMesh, 
        std::vector<Vector3D> &holeCenters,
        tetgenio &out
    ){
        out.numberofpoints = surfaceMesh.nodes.size();
        out.firstnumber = 0;

        out.pointlist = new REAL[out.numberofpoints*3];
        out.pointmarkerlist = new int[out.numberofpoints];
        for(int i=0; i<out.numberofpoints; i++){
            out.pointlist[3*i] = surfaceMesh.nodes[i]->pos[0];
            out.pointlist[3*i+1] = surfaceMesh.nodes[i]->pos[1];
            out.pointlist[3*i+2] = surfaceMesh.nodes[i]->pos[2];
            out.pointmarkerlist[i] = surfaceMesh.nodes[i]->label;
        }

        out.numberoffacets = surfaceMesh.triangles.size();
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
            p->vertexlist[0] = surfaceMesh.triangles[i]->nodes[0]->index;
            p->vertexlist[1] = surfaceMesh.triangles[i]->nodes[1]->index;
            p->vertexlist[2] = surfaceMesh.triangles[i]->nodes[2]->index;

            out.facetmarkerlist[i] = surfaceMesh.triangles[i]->label;
        }
    	out.numberofholes = holeCenters.size();
	    out.holelist = new double[3*out.numberofholes];
	    for(int i=0; i<holeCenters.size(); i++){
            out.holelist[3*i    ] = holeCenters[i][0];
            out.holelist[3*i + 1] = holeCenters[i][1];
            out.holelist[3*i + 2] = holeCenters[i][2];
        }
    }

    void transportTETGENIOToMesh(tetgenio &in, Mesh &out, int withLabel){
        int numNodes = in.numberofpoints;
        int numTets = in.numberoftetrahedra;
        for(int i=0; i<numNodes; i++){
            Node *n = out.addNode(&(in.pointlist[i*3])); 
            if (withLabel){
                n->label = in.pointmarkerlist[i];
            }
            
        }

        int first = in.firstnumber;
        for(int i=0; i<numTets; i++){
            int base = i*4;
            Tetrahedron *tet = out.addTetrahedron(out.nodes[in.tetrahedronlist[base]-first]
            , out.nodes[in.tetrahedronlist[base+1] - first]
            , out.nodes[in.tetrahedronlist[base+2] - first]
            , out.nodes[in.tetrahedronlist[base+3] - first]);
        }

        out.rebuildIndices();        
    }

}