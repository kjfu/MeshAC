#include "mesher3d_core.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include "mesher3d_io.h"
#include <iostream>
#include <ctime>
#include "mesh.h"
#include "triangle.h"
#include "SurfaceMesh.h"
#include "MeshRefiner.h"
#include <unordered_set>
#include <fstream>
#include "TetgenTool.h"
#include "TriangleTool.h"
namespace MeshAC{
void generateConvexHull(const std::string &fileIn, const std::string &fileOut){
	tetgenio in, out;
	loadMesh(&in, fileIn);
	char cmd[] ="Q";
	tetrahedralize(cmd, &in, &out);
	saveAsMESH(&out, fileOut);
}


//dd
void generateBoundingBoxTETGENIO(tetgenio &tetIn, Vector3D xyzmax, Vector3D xyzmin, double size, tetgenio &tetOut){

	tetIn.firstnumber=0;
	tetIn.numberofpoints = 8;
	tetIn.pointlist = new double[3*tetIn.numberofpoints];
	tetIn.pointlist[0]  = xyzmin[0]; tetIn.pointlist[1]  = xyzmin[1]; tetIn.pointlist[2]  = xyzmin[2];
	tetIn.pointlist[3]  = xyzmax[0]; tetIn.pointlist[4]  = xyzmin[1]; tetIn.pointlist[5]  = xyzmin[2];
	tetIn.pointlist[6]  = xyzmax[0]; tetIn.pointlist[7]  = xyzmax[1]; tetIn.pointlist[8]  = xyzmin[2];
	tetIn.pointlist[9]  = xyzmin[0]; tetIn.pointlist[10] = xyzmax[1]; tetIn.pointlist[11] = xyzmin[2];
	tetIn.pointlist[12] = xyzmin[0]; tetIn.pointlist[13] = xyzmin[1]; tetIn.pointlist[14] = xyzmax[2];
	tetIn.pointlist[15] = xyzmax[0]; tetIn.pointlist[16] = xyzmin[1]; tetIn.pointlist[17] = xyzmax[2];
	tetIn.pointlist[18] = xyzmax[0]; tetIn.pointlist[19] = xyzmax[1]; tetIn.pointlist[20] = xyzmax[2];
	tetIn.pointlist[21] = xyzmin[0]; tetIn.pointlist[22] = xyzmax[1]; tetIn.pointlist[23] = xyzmax[2];
	tetIn.numberoffacets = 6;
    tetIn.facetlist = new tetgenio::facet[tetIn.numberoffacets];
	tetIn.numberofpointmtrs=1;
	tetIn.pointmtrlist = new double[8];
	for(int i=0;i<8; i++){
		tetIn.pointmtrlist[i] = size;
	}

	static int faceIndices[6][4]={
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{1, 2, 6, 5},
		{2, 3, 7, 6},
		{0, 3, 7, 4},
		{1, 0, 4, 5}};

    for(int i=0; i<tetIn.numberoffacets; i++){
        tetgenio::facet *f = &tetIn.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = nullptr;
        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = faceIndices[i][0];
        p->vertexlist[1] = faceIndices[i][1];
        p->vertexlist[2] = faceIndices[i][2];
        p->vertexlist[3] = faceIndices[i][3];
    }

	char cmd[] ="pqmDQ";

	tetrahedralize(cmd, &tetIn, &tetOut);


}


void delaunayTetrahedralizationWithHoles(const std::string &fileIn, const std::string &fileOut, double ssize, bool beQuiet){
	Vector3D xyzmax;
	Vector3D xyzmin;
	// Vector3D oxyzmax;
	// Vector3D oxyzmin;

	SurfaceMesh shellSurface;
	Mesh innerMesh;
	tetgenio innerIn;
	tetgenio innerOut;
	loadNodesWithLabel(innerIn, fileIn, xyzmax, xyzmin);
	Vector3D dxyz=xyzmax-xyzmin;
	double size = fmax(dxyz[0], fmax(dxyz[1], dxyz[2]));
	size = ssize<1e-9 ? 0.1*size: ssize;
	// char cmd[] = "Q";
	// tetrahedralize(cmd, &innerIn, &innerOut);
	generateConvaxHullFromPointsIn3DRemoveHoles(innerIn, innerMesh, shellSurface);
	
	//shellSurface.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/two_holes2/surf.vtk");
	shellSurface.estimateSizing();
	int numNodesInsideShell = shellSurface.nodes.size();



	// tetgenio outerIn;
	// tetgenio outerOut;
	// generateBoundingBoxTETGENIO(outerIn, xyzmax, xyzmin, size, outerOut);

	SurfaceMesh tmp;
	// tmp.loadTETGENIO(outerOut);
	generateBoundingBoxTETGENIO(xyzmax, xyzmin, size,tmp);
	int numNodesOutsideShell = tmp.nodes.size();

	shellSurface.addSubSurfaceMesh(tmp);
	tetgenio shellIn;
	tetgenio shellOut;
	shellSurface.exportTETGENIO(shellIn);

	shellIn.numberofpointmtrs = 1;
	shellIn.pointmtrlist = new double[shellIn.numberofpoints];
	for(int i=0; i<numNodesInsideShell; i++){
		shellIn.pointmtrlist[i] = shellSurface.maxSizing;
	}
	for(int i=numNodesInsideShell; i<shellSurface.nodes.size(); i++){
		shellIn.pointmtrlist[i] = size;
	}
	std::vector<Vector3D> holes;
	innerMesh.getSubRegionCenters(holes);
	shellIn.numberofholes=holes.size();
	shellIn.holelist = new double[3*shellIn.numberofholes];
	for(int i=0; i<holes.size(); i++){
		shellIn.holelist[3*i] = holes[i][0];
		shellIn.holelist[3*i+1] = holes[i][1];
		shellIn.holelist[3*i+2] = holes[i][2];
}
	// Vector3D holePos(0, 0, 0);
	// for(auto n: shellSurface.nodes){
	// 	holePos+=n->pos;
	// }
	// holePos/=shellSurface.nodes.size();
	// shellIn.holelist[0] = holePos[0];
	// shellIn.holelist[1] = holePos[1];
	// shellIn.holelist[2] = holePos[2];
	char cmdcmd[] = "pq1.1/10mYQ";
	tetrahedralize(cmdcmd, &shellIn, &shellOut);
	Mesh shellMesh;
	shellMesh.loadTETGENIO(shellOut);
	for(int i=0; i<numNodesInsideShell; i++){
		shellMesh.nodes[i]->label = 2;
	}
	for(int i=numNodesInsideShell; i<numNodesInsideShell+numNodesOutsideShell; i++){
		shellMesh.nodes[i]->label = 1;
	}
	for( int i=numNodesInsideShell+numNodesOutsideShell; i<shellMesh.nodes.size(); i++){
		shellMesh.nodes[i]->label = 3;
	}
	for(auto &tet: shellMesh.tetrahedrons){
		tet->label=1;
	}

	shellMesh.mergeMesh(innerMesh, shellSurface.minSizing*0.01);
	// shellMesh.exportVTK(fileOut+"o.vtk");
	shellMesh.exportMESH(fileOut);
	shellMesh.exportVTK(fileOut+"2vtk.vtk");
}



void delaunayTetrahedralization(const std::string &fileIn, const std::string &fileOut, double ssize, bool beQuiet){
	Vector3D xyzmax;
	Vector3D xyzmin;
	// Vector3D oxyzmax;
	// Vector3D oxyzmin;

	SurfaceMesh shellSurface;
	Mesh innerMesh;
	tetgenio innerIn;
	tetgenio innerOut;
	loadNodesWithLabel(innerIn, fileIn, xyzmax, xyzmin);
	Vector3D dxyz=xyzmax-xyzmin;
	double size = fmax(dxyz[0], fmax(dxyz[1], dxyz[2]));
	size = ssize<1e-9 ? 0.1*size: ssize;
	// char cmd[] = "Q";
	// tetrahedralize(cmd, &innerIn, &innerOut);
	generateConvaxHullFromPointsIn3D(innerIn, innerMesh, shellSurface);
	
	//shellSurface.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/paper_test/refine_pipleline/atomboundary.vtk");
	shellSurface.estimateSizing();
	int numNodesInsideShell = shellSurface.nodes.size();



	// tetgenio outerIn;
	// tetgenio outerOut;
	// generateBoundingBoxTETGENIO(outerIn, xyzmax, xyzmin, size, outerOut);

	SurfaceMesh tmp;
	// tmp.loadTETGENIO(outerOut);
	generateBoundingBoxTETGENIO(xyzmax, xyzmin, size, tmp);
	int numNodesOutsideShell = tmp.nodes.size();
	//tmp.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/paper_test/refine_pipleline/outboundary.vtk");
	shellSurface.addSubSurfaceMesh(tmp);
	tetgenio shellIn;
	tetgenio shellOut;
	shellSurface.exportTETGENIO(shellIn);

	shellIn.numberofpointmtrs = 1;
	shellIn.pointmtrlist = new double[shellIn.numberofpoints];
	for(int i=0; i<numNodesInsideShell; i++){
		shellIn.pointmtrlist[i] = shellSurface.maxSizing;
	}
	for(int i=numNodesInsideShell; i<shellSurface.nodes.size(); i++){
		shellIn.pointmtrlist[i] = size;
	}
	shellIn.numberofholes=1;
	shellIn.holelist = new double[3];
	Vector3D holePos(0, 0, 0);
	for(auto n: shellSurface.nodes){
		holePos+=n->pos;
	}
	holePos/=shellSurface.nodes.size();
	shellIn.holelist[0] = holePos[0];
	shellIn.holelist[1] = holePos[1];
	shellIn.holelist[2] = holePos[2];
	char cmdcmd[] = "pq1.1/10mYQ";
	tetrahedralize(cmdcmd, &shellIn, &shellOut);
	Mesh shellMesh;
	shellMesh.loadTETGENIO(shellOut);
	for(int i=0; i<numNodesInsideShell; i++){
		shellMesh.nodes[i]->label = 2;
	}
	for(int i=numNodesInsideShell; i<numNodesInsideShell+numNodesOutsideShell; i++){
		shellMesh.nodes[i]->label = 1;
	}
	for( int i=numNodesInsideShell+numNodesOutsideShell; i<shellMesh.nodes.size(); i++){
		shellMesh.nodes[i]->label = 3;
	}
	for(auto &tet: shellMesh.tetrahedrons){
		tet->label=1;
	}

	shellMesh.mergeMesh(innerMesh, shellSurface.minSizing*0.01);
	// shellMesh.exportVTK(fileOut+"o.vtk");
	shellMesh.exportMESH(fileOut);
	shellMesh.exportVTK(fileOut+"2vtk.vtk");


}
void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers, bool beQuiet){

	tetgenio tmp_inside_in;
	tetgenio tmp_inside_out;

	tetgenio tmp_outside_in;
	tetgenio tmp_outside_out;

	std::vector<int> insidePointIndices;
	std::vector<int> outsidePointIndices;
	int np = in->numberofpoints;

	for (int i=0; i< in->numberofpoints; i++){
		if (in->pointmarkerlist[i] == 0){
			insidePointIndices.push_back(i);
		}
		else{
			outsidePointIndices.push_back(i);
		}
	}

	/**
	 * @brief copy points  with specific index to another tetgenio
	 * 
	 */
	auto copyPoint = 
	[&in]
	(auto &tetgenContainer, auto &indexContainer, double *center){
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;
		tetgenContainer.numberofpoints = static_cast<int>(indexContainer.size());
		tetgenContainer.pointlist = new REAL[tetgenContainer.numberofpoints*3];
		tetgenContainer.pointmarkerlist = new int[tetgenContainer.numberofpoints];
		for (int i=0; i<tetgenContainer.numberofpoints; i++){
			tetgenContainer.pointlist[i*3+0] = in->pointlist[indexContainer[i]*3+0];
			tetgenContainer.pointlist[i*3+1] = in->pointlist[indexContainer[i]*3+1];
			tetgenContainer.pointlist[i*3+2] = in->pointlist[indexContainer[i]*3+2];
			center[0] += tetgenContainer.pointlist[i*3+0];
			center[1] += tetgenContainer.pointlist[i*3+1];
			center[2] += tetgenContainer.pointlist[i*3+2];
			tetgenContainer.pointmarkerlist[i] = in->pointmarkerlist[indexContainer[i]];
		}
		center[0] /= static_cast<double>(tetgenContainer.numberofpoints);
		center[1] /= static_cast<double>(tetgenContainer.numberofpoints);
		center[2] /= static_cast<double>(tetgenContainer.numberofpoints);
	};

	double center[3];	
	copyPoint(tmp_outside_in, outsidePointIndices, center);
	copyPoint(tmp_inside_in, insidePointIndices, center);


	/**
	 * @brief construct convex hull of points
	 * 
	 */
	tetgenbehavior b0;
	if(beQuiet){
		char command0[] = "fQ";
		b0.parse_commandline(command0);
	}
	else{
		char command0[] = "f";
		b0.parse_commandline(command0);
	}

	fprintf(stdout, "*****Tetrahedralize inside points!\n");
	tetrahedralize(&b0, &tmp_inside_in, &tmp_inside_out);
	fprintf(stdout, "*****Tetrahedralize outside points!\n");
	tetrahedralize(&b0, &tmp_outside_in, &tmp_outside_out);

	/**
	 * @brief Get the max edge in a tethedra
	 * 
	 */
	double h = 0;
	auto maxEdgeSizeInTet = 
	[&tmp_inside_out]
	(auto firstpoint){
		double hmax =0;
		for(int i=0; i<4; i++){
			int p0 = firstpoint[i]*3;
			int p1 = firstpoint[(i+1)%4]*3;
			double *point = tmp_inside_out.pointlist;
			double length = 
			sqrt(pow((point[p0+0] - point[p1+0]),2)
			+ pow((point[p0+1] - point[p1+1]),2)
			+ pow((point[p0+2] - point[p1+2]),2));
			hmax = hmax<length ? length : hmax;
		}
		return hmax;
	}; 

	for (int i=0; i<tmp_inside_out.numberoftetrahedra; i++){
		double hmax = maxEdgeSizeInTet(&tmp_inside_out.tetrahedronlist[4*i]);
		h = h < hmax ? hmax : h; 
	}
	
	size = size > h ? size : h;

	/**
	 * @brief Count boundary triangle facets on a convex hull
	 * 
	 */

	auto travelBoundaryTriangleFacets =
	[]
	(auto &tetgenContainer, auto &pointIndexContainer, auto &facetContainer){				
		std::map<int, int> indexmap;
		for(int i=0; i<tetgenContainer.numberoftrifaces; i++){
			if (tetgenContainer.trifacemarkerlist[i]==1){
				facetContainer.emplace_back();
				for (int j=0; j<3; j++){
					if (std::find(pointIndexContainer.begin(), pointIndexContainer.end(), tetgenContainer.trifacelist[3*i+j]) == pointIndexContainer.end()){
					
						indexmap[tetgenContainer.trifacelist[3*i+j]] = pointIndexContainer.size();
						facetContainer.back().push_back(pointIndexContainer.size());						
						pointIndexContainer.push_back(tetgenContainer.trifacelist[3*i+j]);
						
					}
					else{
						facetContainer.back().push_back(indexmap[tetgenContainer.trifacelist[3*i+j]]);
					}
				}

			}
		}
	};


	tetgenio tmp_shell_in;
	tetgenio tmp_shell_out;
	std::vector<int> outsideShellPointIndices;
	std::vector<int> insideShellPointIndices;
	std::vector<std::vector<int>> outsideFacets;
	std::vector<std::vector<int>> insideFacets;



	travelBoundaryTriangleFacets(tmp_outside_out, outsideShellPointIndices, outsideFacets);

	/**
	 * @brief  Refine outside mesh
	 * 
	 */

	tetgenio tmp_outside_refine_in;
	tetgenio tmp_outside_refine_out;
	tmp_outside_refine_in.numberoffacets = outsideFacets.size();
	tmp_outside_refine_in.facetlist = new tetgenio::facet[tmp_outside_refine_in.numberoffacets];
	tmp_outside_refine_in.numberofholes = 0;
	tmp_outside_refine_in.facetmarkerlist = new int[tmp_outside_refine_in.numberoffacets];
	for(int i=0; i<tmp_outside_refine_in.numberoffacets; i++){
		tmp_outside_refine_in.facetmarkerlist[i]=1;
		tetgenio::facet &f = tmp_outside_refine_in.facetlist[i];		
		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
		f.numberofholes = 0;
		f.holelist = NULL;
		tetgenio::polygon &p = f.polygonlist[0];
		p.numberofvertices = 3;
		p.vertexlist = new int[p.numberofvertices];
		p.vertexlist[0] = outsideFacets[i][0];
		p.vertexlist[1] = outsideFacets[i][1];
		p.vertexlist[2] = outsideFacets[i][2];	

	}
	tmp_outside_refine_in.numberofpoints = outsideShellPointIndices.size();
	tmp_outside_refine_in.pointlist = new REAL[3*tmp_outside_refine_in.numberofpoints];		
	tmp_outside_refine_in.pointmarkerlist = new int[tmp_outside_refine_in.numberofpoints];
	auto index = 0;		
	for(auto i: outsideShellPointIndices){
		tmp_outside_refine_in.pointlist[(index)*3 + 0] = tmp_outside_out.pointlist[3*i + 0];
		tmp_outside_refine_in.pointlist[(index)*3 + 1] = tmp_outside_out.pointlist[3*i + 1];
		tmp_outside_refine_in.pointlist[(index)*3 + 2] = tmp_outside_out.pointlist[3*i + 2];
		tmp_outside_refine_in.pointmarkerlist[index] = tmp_outside_out.pointmarkerlist[i];
		index++;
	}

	tmp_outside_refine_in.numberofpointmtrs = 1;
	tmp_outside_refine_in.pointmtrlist = new REAL[tmp_outside_refine_in.numberofpoints];
	for (int i=0; i<tmp_outside_refine_in.numberofpoints; i++){
		tmp_outside_refine_in.pointmtrlist[i] = size;
	}
	tetgenbehavior b0_refine;
	if(beQuiet){
		char command_refine[] = "pqmDfQ";
		b0_refine.parse_commandline(command_refine);		
	}
	else{
		char command_refine[] = "pqmDf";
		b0_refine.parse_commandline(command_refine);	
	}

	fprintf(stdout, "*****Refine outside!\n");
	tetrahedralize(&b0_refine, &tmp_outside_refine_in, &tmp_outside_refine_out);
	/**
	 * @brief  Build a shell with tetgenio
	 * 
	 */
	outsideShellPointIndices.clear();
	outsideFacets.clear();
	travelBoundaryTriangleFacets(tmp_inside_out, insideShellPointIndices, insideFacets);
	travelBoundaryTriangleFacets(tmp_outside_refine_out, outsideShellPointIndices, outsideFacets);
	
	tmp_shell_in.numberofpoints = insideShellPointIndices.size() + outsideShellPointIndices.size();
	tmp_shell_in.pointlist = new REAL[tmp_shell_in.numberofpoints*3];
	tmp_shell_in.numberofpointmtrs = 1;
	tmp_shell_in.pointmarkerlist = new int[tmp_shell_in.numberofpoints];

	tmp_shell_in.numberoffacets = insideFacets.size() + outsideFacets.size();
	tmp_shell_in.facetlist = new tetgenio::facet[tmp_shell_in.numberoffacets];
	tmp_shell_in.facetmarkerlist = new int[tmp_shell_in.numberoffacets];
	tmp_shell_in.numberofholes = 1;
	tmp_shell_in.holelist = new double[3];
	tmp_shell_in.holelist[0] = center[0];
	tmp_shell_in.holelist[1] = center[1];
	tmp_shell_in.holelist[2] = center[2];

	auto addPoints =
	[&tmp_shell_in]
	(auto &tetgenioContainer, auto &pointIndexContainer, auto base){
		auto index = base;		
		for(auto i: pointIndexContainer){
			tmp_shell_in.pointlist[(index)*3 + 0] = tetgenioContainer.pointlist[3*i + 0];
			tmp_shell_in.pointlist[(index)*3 + 1] = tetgenioContainer.pointlist[3*i + 1];
			tmp_shell_in.pointlist[(index)*3 + 2] = tetgenioContainer.pointlist[3*i + 2];
			tmp_shell_in.pointmarkerlist[index] = tetgenioContainer.pointmarkerlist[i];
			index++;
		}
		// tmp_shell_in.pointmarkerlist[index] = 1;
	};
	addPoints(tmp_inside_out, insideShellPointIndices, 0);
	addPoints(tmp_outside_refine_out, outsideShellPointIndices, insideShellPointIndices.size());



	auto addTriangleFacet = 
	[&tmp_shell_in]
	(auto &tetgenioContainer, auto &facetContainer, int facetBase, int pointBase){
		for(int i=0; i<facetContainer.size(); i++){
			tetgenio::facet &f = tmp_shell_in.facetlist[i+facetBase];		
			f.numberofpolygons = 1;
			f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
			f.numberofholes = 0;
			f.holelist = NULL;
			tetgenio::polygon &p = f.polygonlist[0];
			p.numberofvertices = 3;
			p.vertexlist = new int[p.numberofvertices];
			p.vertexlist[0] = facetContainer[i][0] + pointBase;
			p.vertexlist[1] = facetContainer[i][1] + pointBase;
			p.vertexlist[2] = facetContainer[i][2] + pointBase;
		}
	};
	addTriangleFacet(tmp_inside_out, insideFacets, 0, 0);
	addTriangleFacet(tmp_outside_refine_out, outsideFacets, insideFacets.size(), insideShellPointIndices.size());


	tmp_shell_in.numberofpointmtrs = 1;
	tmp_shell_in.pointmtrlist = new REAL[tmp_shell_in.numberofpoints];
    for(int i=0; i<insideShellPointIndices.size(); i++){

        tmp_shell_in.pointmtrlist[i] = h;
    }
	for(int i=0; i<outsideShellPointIndices.size(); i++){
		tmp_shell_in.pointmtrlist[insideShellPointIndices.size()+i] = size;
	}
	
	tetgenbehavior b1;
	if (beQuiet){
		char command1[] = "pqmYQ";
		b1.parse_commandline(command1);
	}
	else{
		char command1[] = "pqmY";
		b1.parse_commandline(command1);
	}

	fprintf(stdout, "*****Tetrahedralize shell!\n");
	tetrahedralize(&b1, &tmp_shell_in, &tmp_shell_out);
	out->firstnumber = 1;
	out->numberoftetrahedra = tmp_shell_out.numberoftetrahedra + tmp_inside_out.numberoftetrahedra;
	out->tetrahedronlist =  new int[out->numberoftetrahedra*4];
	out->numberofpoints = tmp_shell_out.numberofpoints + tmp_inside_out.numberofpoints - insideShellPointIndices.size();

	out->pointlist = new REAL[out->numberofpoints * 3];
	out->pointmarkerlist = new int[out->numberofpoints];	
	for(int i=0; i<tmp_inside_out.numberofpoints; i++){
		out->pointlist[3*i + 0] = tmp_inside_out.pointlist[3*i + 0];
		out->pointlist[3*i + 1] = tmp_inside_out.pointlist[3*i + 1];
		out->pointlist[3*i + 2] = tmp_inside_out.pointlist[3*i + 2];
		out->pointmarkerlist[i] = 0;
	}	

	int minusPointBase = 0 + tmp_inside_out.numberofpoints - insideShellPointIndices.size(); 
	for(int i=insideShellPointIndices.size(); i<tmp_shell_out.numberofpoints; i++){
	
		out->pointlist[3*(i+minusPointBase) + 0] = tmp_shell_out.pointlist[3*i + 0];
		out->pointlist[3*(i+minusPointBase) + 1] = tmp_shell_out.pointlist[3*i + 1];
		out->pointlist[3*(i+minusPointBase) + 2] = tmp_shell_out.pointlist[3*i + 2];
		out->pointmarkerlist[i+minusPointBase] = 1;
	}

	for(int i=0; i<tmp_inside_out.numberoftetrahedra; i++){
		out->tetrahedronlist[4*i + 0] = tmp_inside_out.tetrahedronlist[4*i + 0]+1;
		out->tetrahedronlist[4*i + 1] = tmp_inside_out.tetrahedronlist[4*i + 1]+1;
		out->tetrahedronlist[4*i + 2] = tmp_inside_out.tetrahedronlist[4*i + 2]+1;
		out->tetrahedronlist[4*i + 3] = tmp_inside_out.tetrahedronlist[4*i + 3]+1;
		tetMarkers.push_back(0);
	}
	
	tetgenio final_convex_out;
	tetgenbehavior b2;
	if (beQuiet){
		char command2[] = "fQ";
		b2.parse_commandline(command2);
	}
	else{
		char command2[] = "f";
		b2.parse_commandline(command2);		
	}

	//printf("asdadasdasdads %d\n", tmp_shell_out.numberofpoints);
	fprintf(stdout, "*****Tetrahedralize all points!\n");
	tetrahedralize(&b2, &tmp_shell_out, &final_convex_out);
	for(int i=0; i<final_convex_out.numberoftrifaces; i++){
		if (final_convex_out.trifacemarkerlist[i]==1){
			for(int j=0; j<3; j++){
				out->pointmarkerlist[final_convex_out.trifacelist[3*i+j]+minusPointBase] = 2;
			}
		}
	}



	auto updateIndex=
	[minusPointBase, &insideShellPointIndices]
	(auto pointIndex){
		if (pointIndex<insideShellPointIndices.size()){
			return insideShellPointIndices[pointIndex] + 1;
		}
		else{
			return pointIndex+minusPointBase+1;
		}
	};



	for (int i=0; i<tmp_shell_out.numberoftetrahedra; i++){
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 0]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+0]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 1]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+1]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 2]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+2]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 3]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+3]);
		tetMarkers.push_back(1);
	}	
	
	// int sumsum=0;
	// int queryIndex= 1;
	// for (int i=0; i<out->numberoftetrahedra; i++){
	// 	if (out->tetrahedronlist[4*i+0]==queryIndex || out->tetrahedronlist[4*i+1]==queryIndex  || out->tetrahedronlist[4*i+2]==queryIndex || out->tetrahedronlist[4*i+3]==queryIndex){
	// 		fprintf(stdout, "***%d  %d  %d  %d\n", out->tetrahedronlist[4*i+0], out->tetrahedronlist[4*i+1], out->tetrahedronlist[4*i+2], out->tetrahedronlist[4*i+3]);
	// 		sumsum++;
	// 	}
	// }
	// fprintf(stdout, "********************* %d\n", sumsum);
}






void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size, bool beQuiet){
	in->numberofpointmtrs = 1;
    in->pointmtrlist = new REAL[in->numberofpoints];
    for(int i=0; i<in->numberofpoints; i++){

        in->pointmtrlist[i] = size;
    }

	tetgenbehavior b;
	if(beQuiet){
		char command[] = "pqmDQ";
		b.parse_commandline(command);		
	}
	else{
		char command[] = "pqmD";
		b.parse_commandline(command);			
	}

    tetrahedralize(&b, in, out);

}

void extractBorder(std::vector<Tetrahedron *>&tets, SurfaceMesh &aSurface){
    std::set<Node *> nodeSet; 
    // HashFacetTable facetTable;
	std::unordered_set<SubTriangle, SubTriangleHasher, SubTriangleEqual> subTriangleSet;
    for(auto e: tets){
        for(int i=0; i<4; i++){
            SubTriangle keyFacet = e->getSubTriangle(i);
            // TriangleFacet goalFacet;
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

		Triangle *tri = new Triangle(getNode(f.forms[0]), getNode(f.forms[1]), getNode(f.forms[2]));
		aSurface.triangles.push_back(tri);
        
    }
	for(auto kv: oldNewNodes){
		aSurface.nodes.push_back(kv.second);
	}

	aSurface.rebuildIndices();
}

void refineMeshV3(const std::string &fileInHead, const std::string &fileOutHead, double hmax, double hmin){
	Mesh goalMesh;
	Mesh backgroundMesh;
	double timebegin, timeend;
	std::vector<int> refine_elements;
	std::vector<std::array<double,3>> append_points;
	goalMesh.loadMESH(fileInHead + ".mesh");
	loadREMESH(refine_elements, append_points, fileInHead+".remesh");
	backgroundMesh.clone(goalMesh);
	backgroundMesh.loadNodeValues(fileInHead + ".value");
	timebegin= clock();
	// std::vector<Vector3D> positions;
	// 	for(auto nt:refine_elements){
	// 		Vector3D vec = goalMesh.tetrahedrons[nt-1]->center();
	// 		positions.push_back(vec);
	// 	}
	// timebegin = clock_t();
	// goalMesh.CavityBasedInsert(positions);
	// timeend = clock_t();
	// std::cout << "[*************] Insert time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
	// goalMesh.exportVTK(fileOutHead + ".1.vtk");

	std::vector<Vector3D> points;
	for(auto &n: goalMesh.nodes){
		if(n->label==2 || n->label==0){
			points.emplace_back(n->pos);
		}
	}
	for(auto &p: append_points){
		points.emplace_back(p[0], p[1], p[2]);
	}

	tetgenio tmpOut;
	tetgenio tmpIn;
	transportVector3dsToTETGENIO(points, tmpIn);
	Mesh innerMesh;
	SurfaceMesh surfaceInside;
	generateConvaxHullFromPointsIn3DRemoveHoles(tmpIn, innerMesh, surfaceInside);

	innerMesh.readyForSpatialSearch();
	// goalMesh.readyForSpatialSearch();






	SurfaceMesh surfaceOutside;
	goalMesh.extractBorder(surfaceOutside);
	surfaceOutside.estimateSizing();
	surfaceInside.estimateSizing();
	int numOutside = surfaceOutside.nodes.size();
	surfaceOutside.addSubSurfaceMesh(surfaceInside);
	int numInside = surfaceInside.nodes.size();

	tetgenio shellIn;
	tetgenio shellOut;
	Mesh shell;
	surfaceOutside.exportTETGENIO(shellIn);

	shellIn.numberofpointmtrs = 1;
	shellIn.pointmtrlist = new double[shellIn.numberofpoints];
	for(int i=0; i<numOutside; i++){
		shellIn.pointmtrlist[i] = hmax>0? hmax:surfaceOutside.minSizing;
	}
	for(int i=numOutside; i<numInside+numOutside; i++){
		shellIn.pointmtrlist[i] = hmin>0? hmin:surfaceInside.minSizing;
	}


	shellIn.numberofholes = 1;
	Vector3D hole;
	for(auto n: surfaceInside.nodes){
		hole+=n->pos;
	}
	hole/=surfaceInside.nodes.size();
	shellIn.holelist = new double[3];
	shellIn.holelist[0] = hole[0];
	shellIn.holelist[1] = hole[1];
	shellIn.holelist[2] = hole[2];
	char cmd[] = "pqmYQ";
	tetrahedralize(cmd, &shellIn, &shellOut);
	Mesh midMesh;
	midMesh.loadTETGENIO(shellOut, true);
	for(auto &tet: innerMesh.tetrahedrons){
		tet->label = 0;
	}
	for(auto &tet: midMesh.tetrahedrons){
		tet->label = 1;
	}

	
	for(int i=surfaceOutside.nodes.size(); i<midMesh.nodes.size(); i++){
		midMesh.nodes[i]->label = 3;
	}

	innerMesh.mergeMesh(midMesh, 1e-10);

	timeend = clock();
	std::cout <<"[*************] Remesh time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
	

	timebegin= clock();
	backgroundMesh.interpolateNodeValuesForAnotherMesh(innerMesh);
	timeend = clock();
	std::cout <<"[*************] Interpolation time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;

	innerMesh.exportNodeValues(fileOutHead + ".value");
	innerMesh.exportMESH(fileOutHead + ".mesh");
	innerMesh.exportVTK(fileOutHead+".vtk");
	std::cout << "Finish Adaption!\n";	
}

void refineMeshV2(const std::string &fileInHead, const std::string &fileOutHead, bool beQuiet){
	Mesh goalMesh;
	Mesh backgroundMesh;
	double timebegin, timeend;
	std::vector<int> refine_elements;
	std::vector<std::array<double,3>> append_points;
	goalMesh.loadMESH(fileInHead + ".mesh");
	loadREMESH(refine_elements, append_points, fileInHead+".remesh");
	backgroundMesh.clone(goalMesh);
	backgroundMesh.loadNodeValues(fileInHead + ".value");
	// std::vector<Vector3D> positions;
	// 	for(auto nt:refine_elements){
	// 		Vector3D vec = goalMesh.tetrahedrons[nt-1]->center();
	// 		positions.push_back(vec);
	// 	}
	MeshRefiner aRefiner(&goalMesh);
	timebegin = clock_t();
	//goalMesh.CavityBasedInsert(positions);
	aRefiner.refine(1,refine_elements, 1, 3);
	timeend = clock_t();
	std::cout << "[*************] Insert time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
	// goalMesh.exportVTK(fileOutHead + ".1.vtk");

	std::vector<Vector3D> points;
	for(auto &n: goalMesh.nodes){
		if(n->label==2 || n->label==0){
			points.emplace_back(n->pos);
		}
	}
	for(auto &p: append_points){
		points.emplace_back(p[0], p[1], p[2]);
	}

	tetgenio tmpOut;
	tetgenio tmpIn;
	transportVector3dsToTETGENIO(points, tmpIn);
	Mesh innerMesh;
	SurfaceMesh surfaceInside;
	generateConvaxHullFromPointsIn3DRemoveHoles(tmpIn, innerMesh, surfaceInside);

	innerMesh.readyForSpatialSearch();
	// goalMesh.readyForSpatialSearch();

	timebegin= clock();
	goalMesh.checkBooleanRemoveSpecial(innerMesh, 0);



	for(auto n: goalMesh.nodes){
		n->edit=1;
	}
	std::vector<Tetrahedron *> rmtets;
	SurfaceMesh surfaceOutside;
	for(auto &e: goalMesh.tetrahedrons){
		if (e->edit!=-1){

			for(auto &n: e->nodes){
				n->edit=0;
			}
		}
		else{
			rmtets.push_back(e);
		}
	}

	extractBorder(rmtets, surfaceOutside);

	for(int i=0; i<goalMesh.nodes.size(); i++){
		if(goalMesh.nodes[i]->edit==1){
			delete goalMesh.nodes[i];
			goalMesh.nodes[i] = goalMesh.nodes.back();
			goalMesh.nodes.pop_back();
			i--;
		}
	}

	for(int i=0; i<goalMesh.tetrahedrons.size(); i++){
		if(goalMesh.tetrahedrons[i]->edit==-1){
			delete goalMesh.tetrahedrons[i];
			goalMesh.tetrahedrons[i] = goalMesh.tetrahedrons.back();
			goalMesh.tetrahedrons.pop_back();
			i--;
		}
	}

	goalMesh.rebuildIndices();
	timeend = clock();
	std::cout <<"[*************] Boolean remove time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;

	// for(auto n: surfaceInside.nodes){
	// 	n->label = 1;
	// }
	// surfaceOutside.estimateSizing();
	// surfaceInside.estimateSizing();

	// int numInside = surfaceInside.nodes.size();
	// int numOutside = surfaceOutside.nodes.size();

	surfaceOutside.addSubSurfaceMesh(surfaceInside);
	// surfaceOutside.exportVTK(fileOutHead + ".2.vtk");

	tetgenio shellIn;
	tetgenio shellOut;
	Mesh shell;
	surfaceOutside.exportTETGENIO(shellIn);

	// shellIn.numberofpointmtrs = 1;
	// shellIn.pointmtrlist = new double[shellIn.numberofpoints];
	// for(int i=0; i<numOutside; i++){
	// 	shellIn.pointmtrlist[i] = surfaceOutside.maxSizing;
	// }
	// for(int i=numOutside; i<numInside+numOutside; i++){
	// 	shellIn.pointmtrlist[i] = surfaceInside.maxSizing;
	// }


	shellIn.numberofholes = 1;
	Vector3D hole;
	for(auto n: surfaceInside.nodes){
		hole+=n->pos;
	}
	hole/=surfaceInside.nodes.size();
	shellIn.holelist = new double[3];
	shellIn.holelist[0] = hole[0];
	shellIn.holelist[1] = hole[1];
	shellIn.holelist[2] = hole[2];
	char cmd[] = "pq1.1/10mYQ";
	tetrahedralize(cmd, &shellIn, &shellOut);
	Mesh midMesh;
	midMesh.loadTETGENIO(shellOut, true);
	for(auto &tet: innerMesh.tetrahedrons){
		tet->label = 0;
	}
	for(auto &tet: midMesh.tetrahedrons){
		tet->label = 1;
	}
	for(auto &tet: goalMesh.tetrahedrons){
		tet->label = 1;
	}
	
	for(int i=surfaceOutside.nodes.size(); i<midMesh.nodes.size(); i++){
		midMesh.nodes[i]->label = 3;
	}

	innerMesh.mergeMesh(midMesh, 1e-10);


	goalMesh.mergeMesh(innerMesh, 1e-10);
	

	timebegin= clock();
	backgroundMesh.interpolateNodeValuesForAnotherMesh(goalMesh);
	timeend = clock();
	std::cout <<"[*************] Interpolation time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;

	goalMesh.exportNodeValues(fileOutHead + ".value");
	goalMesh.exportMESH(fileOutHead + ".mesh");
	goalMesh.exportVTK(fileOutHead+".vtk");
	std::cout << "Finish Adaption!\n";

}


void refineMesh(const std::string &fileInHead, const std::string &fileOutHead, bool beQuiet){

		Mesh goalMesh;
		Mesh backgroundMesh;
		double timebegin, timeend;
		std::vector<int> refine_elements;
		std::vector<std::array<double,3>> append_points;
		goalMesh.loadMESH(fileInHead+".mesh");		
		loadREMESH(refine_elements, append_points, fileInHead+".remesh");		
		backgroundMesh.clone(goalMesh);
		backgroundMesh.loadNodeValues(fileInHead+".value");	
		// std::cout << "[*************] Input nodes: "<<goalMesh.nodes.size() << "; input tets: "<< goalMesh.tetrahedrons.size() << std::endl;
		std::vector<Vector3D> positions;
		for(auto nt:refine_elements){
			Vector3D vec = goalMesh.tetrahedrons[nt-1]->center();
			positions.push_back(vec);
		}
		timebegin = clock_t();
		goalMesh.CavityBasedInsert(positions);
		
		timeend = clock_t();
		std::cout << "[*************] Insert time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
		// std::cout << "[*************] After insert, nodes: "<<goalMesh.nodes.size() << "; tets: "<< goalMesh.tetrahedrons.size() << std::endl;		






		goalMesh.estimateSizing();
		std::vector<Tetrahedron *> tetsToRemove;
		std::vector<Node *> sNodes;
		for (auto &e: goalMesh.tetrahedrons){
			bool removeOut = true;
			for(auto n: e->nodes){
				if (n->label!=0){
					removeOut = false;
					break;
				}

			}
			if (removeOut){
				tetsToRemove.push_back(e);
			}
		}

		Mesh tmpMesh;
		for (auto n: goalMesh.nodes){
			if (n->label==0){
				Node *node = new Node(*n);
				tmpMesh.nodes.push_back(node);
			}
		}

		for(auto p:append_points){
			Node *node = new Node(p.data());
			tmpMesh.nodes.push_back(node);
		}




		tetgenio tmpout;
		tetgenio tmpin;
		Mesh gradMesh;
		std::vector<Node *> grad_nodes;
		std::vector<SubTriangle> grad_facets; 


		transportNodesToTETGENIO(tmpMesh.nodes, tmpin);
		if(beQuiet){
			char cmd[]="fQ";
			tetrahedralize(cmd, &tmpin, &tmpout);			
		}
		else{
			char cmd[]="f";
			tetrahedralize(cmd, &tmpin, &tmpout);			
		}

		instructTetrahedronConnectByTETGENIO(tmpMesh.nodes, tmpout, tmpMesh.tetrahedrons);		
		tmpMesh.rebuildIndices();
		tmpMesh.estimateSizing();
		tmpMesh.extractBorder(grad_nodes, grad_facets);
		Vector3D hole;
		for(auto n: grad_nodes){
			hole+=n->pos;
		}
		hole/=grad_nodes.size();



		tmpMesh.readyForSpatialSearch();
		timebegin= clock();
		goalMesh.checkBooleanRemove(tmpMesh);
		timeend = clock();


		std::cout <<"[*************] Boolean remove time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;


		std::vector<Tetrahedron *> rmTets;
		for(auto e: goalMesh.tetrahedrons){
			if(e->edit==-1){
				rmTets.push_back(e);
			}
		}
		extratctBorder(rmTets, grad_nodes, grad_facets);


		goalMesh.deleteLargeScaleTetrahedronsPermanently(rmTets);

		std::vector<Vector3D> holes;
		holes.emplace_back(0,0,0);
		tetgenio grad_in, grad_out;
		transportFacetsToTETGENIO(grad_nodes, grad_facets, holes, grad_in);
		if (beQuiet){
			char cmd0[]="pqmYQ";
			tetrahedralize(cmd0, &grad_in, &grad_out);
		}
		else{
			char cmd0[]="pqmY";
			tetrahedralize(cmd0, &grad_in, &grad_out);	
		}

		gradMesh.loadTETGENIO(grad_out);

		std::vector<Node *>mergeNodes0;
		// std::cout << "[*************] Start merge!\n";
		goalMesh.mergeMesh(gradMesh, mergeNodes0);	
		for(auto n:goalMesh.nodes){
			n->label = 1;
		}
		for(auto e:goalMesh.tetrahedrons){
			e->label = 1;
		}


		for(auto n: tmpMesh.nodes){
			n->label = 0;
		}
		for(auto e:tmpMesh.tetrahedrons){
			e->label = 0;
		}
		std::vector<Node *>mergeNodes1;
		goalMesh.mergeMesh(tmpMesh, mergeNodes1);

		for(auto n: mergeNodes1){
			n->label = 0;
		}
		std::vector<Node *> nodes;
		extractBorderNodes(goalMesh.tetrahedrons, nodes);
		for(auto n: nodes){
			n->label = 2;
		}	
		// std::cout << "[*************] Finish merge!\n";
		timebegin= clock();
		backgroundMesh.interpolateNodeValuesForAnotherMesh(goalMesh);
		timeend = clock();
		std::cout <<"[*************] Interpolation time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;

		goalMesh.exportNodeValues(fileOutHead + ".value");
		goalMesh.exportMESH(fileOutHead + ".mesh");
		std::cout << "Finish Adaption!\n";
}

void resetPoints(tetgenio &tet, Vector3D pMax, Vector3D pMin, std::vector<int> &indexOf1){
	//0
	tet.pointlist[3*indexOf1[0]] = pMax[0];
	tet.pointlist[3*indexOf1[0]+1] = pMax[1];
	tet.pointlist[3*indexOf1[0]+2] = pMax[2];
	
	//1
	tet.pointlist[3*indexOf1[1]] = pMin[0];
	tet.pointlist[3*indexOf1[1]+1] = pMax[1];
	tet.pointlist[3*indexOf1[1]+2] = pMax[2];

	//2
	tet.pointlist[3*indexOf1[2]] = pMax[0];
	tet.pointlist[3*indexOf1[2]+1] = pMin[1];
	tet.pointlist[3*indexOf1[2]+2] = pMax[2];

	//3	
	tet.pointlist[3*indexOf1[3]] = pMax[0];
	tet.pointlist[3*indexOf1[3]+1] = pMax[1];
	tet.pointlist[3*indexOf1[3]+2] = pMin[2];

	//4
	tet.pointlist[3*indexOf1[4]] = pMin[0];
	tet.pointlist[3*indexOf1[4]+1] = pMin[1];
	tet.pointlist[3*indexOf1[4]+2] = pMax[2];

	//5
	tet.pointlist[3*indexOf1[5]] = pMax[0];
	tet.pointlist[3*indexOf1[5]+1] = pMin[1];
	tet.pointlist[3*indexOf1[5]+2] = pMin[2];

	//6
	tet.pointlist[3*indexOf1[6]] = pMin[0];
	tet.pointlist[3*indexOf1[6]+1] = pMax[1];
	tet.pointlist[3*indexOf1[6]+2] = pMin[2];

	//7
	tet.pointlist[3*indexOf1[7]] = pMin[0];
	tet.pointlist[3*indexOf1[7]+1] = pMin[1];
	tet.pointlist[3*indexOf1[7]+2] = pMin[2];

}
void generateBoundingBoxTETGENIO(Vector3D xyzmax, Vector3D xyzmin, double size, SurfaceMesh &aSurface){
	std::vector<std::array<int, 2>> bottomEdges;
	std::vector<std::array<double, 2>> bottomEdgeNodes;

	SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
	std::vector<std::array<double,2>> holes;
	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triTopBottom;
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triTopBottom);
	faceTop.projectTRIANGULATEIO(triTopBottom, PROJECTION_TYPE::XY_PLANE, xyzmax[2]);
	faceBottom.projectTRIANGULATEIO(triTopBottom, PROJECTION_TYPE::XY_PLANE, xyzmin[2]);
	deleteTRIANGULATEIOAllocatedArrays(triTopBottom);

	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[1], xyzmax[2]}, {xyzmin[1],xyzmin[2]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triFrontBack;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triFrontBack);
	faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmax[0]);
	faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmin[0]);
	deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[2], xyzmax[0]}, {xyzmin[2],xyzmin[0]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triLeftRight;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triLeftRight);
	faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmax[1]);
	faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmin[1]);
	deleteTRIANGULATEIOAllocatedArrays(triLeftRight);

	double eps = 0.01*size;
	aSurface.mergeSurfaceMesh(faceRight,eps);
	aSurface.mergeSurfaceMesh(faceLeft, eps);
	aSurface.mergeSurfaceMesh(faceFront, eps);
	aSurface.mergeSurfaceMesh(faceBack, eps);
	aSurface.mergeSurfaceMesh(faceTop, eps);	
	aSurface.mergeSurfaceMesh(faceBottom, eps);

}


void parseZHandleV3(SurfaceMesh &zHandleSurface, Vector3D xyzmax, Vector3D xyzmin, Vector3D oxyzmax, Vector3D oxyzmin, double size, Mesh &meshOut){

	std::vector<std::array<int, 2>> bottomEdges;
	std::vector<std::array<double, 2>> bottomEdgeNodes;
	std::vector<std::array<int, 2>> topEdges;
	std::vector<std::array<double, 2>> topEdgeNodes;
	parseZHandle(zHandleSurface, oxyzmax[2], oxyzmin[2], topEdgeNodes, bottomEdgeNodes, topEdges, bottomEdges);

	zHandleSurface.estimateSizing();
	if (size==-1){
		size = zHandleSurface.maxSizing * 10;
	}

	//Generate outer surface mesh	
	SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
	std::array<double, 2> oxymax({oxyzmax[0], oxyzmax[1]});
	std::array<double, 2> oxymin({oxyzmin[0], oxyzmin[1]});

	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, bottomEdgeNodes, bottomEdges);
	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, topEdgeNodes, topEdges);
	std::vector<std::array<double,2>> holes;
	holes.push_back({0.5*(oxymax[0]+oxymin[0]), 0.5*(oxymax[1]+oxymin[1])});
	triangulateio triBottom;
	triangulateio triTop;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triBottom);
	generateMeshInPlaneWithEdges(topEdgeNodes, topEdges, holes, size*size/2, triTop);

	faceBottom.projectTRIANGULATEIO(triBottom, PROJECTION_TYPE::XY_PLANE, xyzmin[2]);
	faceTop.projectTRIANGULATEIO(triTop, PROJECTION_TYPE::XY_PLANE, xyzmax[2]);
	deleteTRIANGULATEIOAllocatedArrays(triBottom);
	deleteTRIANGULATEIOAllocatedArrays(triTop);

	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[1], xyzmax[2]}, {xyzmin[1],xyzmin[2]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triFrontBack;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triFrontBack);
	faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmax[0]);
	faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmin[0]);
	deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[2], xyzmax[0]}, {xyzmin[2],xyzmin[0]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triLeftRight;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triLeftRight);
	faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmax[1]);
	faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmin[1]);
	deleteTRIANGULATEIOAllocatedArrays(triLeftRight);


	auto addNodeLabel=[](auto &mesh, int label){
		for(auto n: mesh.nodes){
			n->label = label;
		}
	};
	addNodeLabel(zHandleSurface, 2);
	addNodeLabel(faceTop, 1);
	addNodeLabel(faceBottom, 1);
	addNodeLabel(faceFront, 1);	
	addNodeLabel(faceBack, 1);	
	addNodeLabel(faceRight, 1);	
	addNodeLabel(faceLeft, 1);
	int numNodes1 = zHandleSurface.nodes.size();
	zHandleSurface.mergeSurfaceMesh(faceRight,1e-10);
	zHandleSurface.mergeSurfaceMesh(faceLeft,1e-10);
	zHandleSurface.mergeSurfaceMesh(faceFront,1e-10);
	zHandleSurface.mergeSurfaceMesh(faceBack,1e-10);
	zHandleSurface.mergeSurfaceMesh(faceTop,1e-10);	
	zHandleSurface.mergeSurfaceMesh(faceBottom,1e-10);
	int numNodes2 = zHandleSurface.nodes.size();

	//Generate tetgenio
	tetgenio in, out;
	zHandleSurface.exportTETGENIO(in);
	in.numberofholes = 1;
	in.holelist = new double[3];

	Vector3D handleHole = Vector3D((xyzmax[0]+xyzmin[0])/2.0, (xyzmax[1]+xyzmin[1])/2.0, (xyzmin[2]+xyzmax[2])/2.0);
	in.holelist[0] = handleHole[0];
	in.holelist[1] = handleHole[1];
	in.holelist[2] = handleHole[2];

	char cmd[] = "pYq1.1/10fQ";
	tetrahedralize(cmd, &in, &out);
	meshOut.loadTETGENIO(out, true);
	//int numTriFacets = out.numberoftrifaces;

	for (int i=0; i<numNodes1; i++){
		meshOut.nodes[i]->label = 2;
	}
	for (int i=numNodes1; i<numNodes2; i++){
		meshOut.nodes[i]->label = 1;
	}
	for (int i=numNodes2; i<meshOut.nodes.size(); i++){
		meshOut.nodes[i]->label = 3;
	}
	for(auto &t: meshOut.tetrahedrons){
		t->label = 1;
	}
}


void parseZHandleV2(SurfaceMesh &zHandleSurface, Vector3D xyzmax, Vector3D xyzmin, Vector3D oxyzmax, Vector3D oxyzmin, double size, Mesh &meshOut){

    std::vector<std::array<Node *, 2>> topEdges;
    std::vector<std::array<Node *, 2>> bottomEdges;
	
	zHandleSurface.rebuildTriangleAdjacency();


	for (auto &n: zHandleSurface.nodes){
		if (abs(n->pos[2]-oxyzmax[2])<1e-13){
			n->edit = 1;
		}
		else if(abs(n->pos[2]-oxyzmin[2])<1e-13){
			n->edit = 2;
		}
		else{
			n->edit = 0;
		}
	}
	std::vector<Triangle *> delTriangles;
	for (auto &t: zHandleSurface.triangles){
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
	for(auto &t: delTriangles){
		for(int i=0; i<3; i++){

			Triangle *tt = t->adjacentTriangles[i];
			if (tt->edit==0){

				if (t->edit==1){
					topEdges.push_back({t->nodes[(i+1)%3],  t->nodes[(i+2)%3]});
				}
				else{
					bottomEdges.push_back({t->nodes[(i+1)%3], t->nodes[(i+2)%3]});						
				}
			}
		}			
	}
	std::vector<Node *> topNodesList;
	std::vector<Node *> bottomNodesList;


	auto clearNodes = []
	(std::vector<std::array<Node *, 2>> &edges, std::vector<Node*> &nodesList){
		std::vector<bool> checked(edges.size(), false);
		while (nodesList.size() != edges.size()){
			if (nodesList.empty()){
				nodesList.push_back(edges[0][1]);
			}
			checked[0] = true; 
			for(int i=1; i<edges.size(); i++){
				if(!checked[i]){
					if(edges[i][0] == nodesList.back()){
						nodesList.push_back(edges[i][1]);
						checked[i] = true;					
					}
					else if(edges[i][1] == nodesList.back()){
						nodesList.push_back(edges[i][0]);
						checked[i] = false;
					}					
				}


			}
		}
	};


	clearNodes(topEdges, topNodesList);
	clearNodes(bottomEdges, bottomNodesList);

	zHandleSurface.deleteTriangles(delTriangles);

	
	zHandleSurface.estimateSizing();
	if (size==-1){
		size = zHandleSurface.maxSizing * 10;
	}
	//Generate tetgenio
	tetgenio in, out;
	int baseNodesIndex = zHandleSurface.nodes.size()*3;
	int numInnerNodes = zHandleSurface.nodes.size();
	int numInnerFacets = zHandleSurface.triangles.size();
	in.numberoffacets = zHandleSurface.triangles.size() + 6;
	in.firstnumber = 0;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];
	in.numberofpoints = zHandleSurface.nodes.size() + 8;
	in.pointlist = new double[3*in.numberofpoints];
	in.pointmtrlist = new double[in.numberofpoints];
	in.pointmarkerlist = new int[in.numberofpoints];
	in.numberofpointmtrs = 1;
	in.numberofholes = 1;
	in.holelist = new double[3];
	for (int i=0; i<zHandleSurface.nodes.size(); i++){
		int base = i*3;
		in.pointlist[base] = zHandleSurface.nodes[i]->pos[0];
		in.pointlist[base+1] = zHandleSurface.nodes[i]->pos[1];
		in.pointlist[base+2] = zHandleSurface.nodes[i]->pos[2];
		in.pointmtrlist[i] = zHandleSurface.maxSizing;
		in.pointmarkerlist[i] = 2;
	}

	in.pointlist[baseNodesIndex+0]  = xyzmin[0]; in.pointlist[baseNodesIndex+1]  = xyzmin[1]; in.pointlist[baseNodesIndex+2]  = xyzmin[2];
	in.pointlist[baseNodesIndex+3]  = xyzmax[0]; in.pointlist[baseNodesIndex+4]  = xyzmin[1]; in.pointlist[baseNodesIndex+5]  = xyzmin[2];
	in.pointlist[baseNodesIndex+6]  = xyzmax[0]; in.pointlist[baseNodesIndex+7]  = xyzmax[1]; in.pointlist[baseNodesIndex+8]  = xyzmin[2];
	in.pointlist[baseNodesIndex+9]  = xyzmin[0]; in.pointlist[baseNodesIndex+10] = xyzmax[1]; in.pointlist[baseNodesIndex+11] = xyzmin[2];
	in.pointlist[baseNodesIndex+12] = xyzmin[0]; in.pointlist[baseNodesIndex+13] = xyzmin[1]; in.pointlist[baseNodesIndex+14] = xyzmax[2];
	in.pointlist[baseNodesIndex+15] = xyzmax[0]; in.pointlist[baseNodesIndex+16] = xyzmin[1]; in.pointlist[baseNodesIndex+17] = xyzmax[2];
	in.pointlist[baseNodesIndex+18] = xyzmax[0]; in.pointlist[baseNodesIndex+19] = xyzmax[1]; in.pointlist[baseNodesIndex+20] = xyzmax[2];
	in.pointlist[baseNodesIndex+21] = xyzmin[0]; in.pointlist[baseNodesIndex+22] = xyzmax[1]; in.pointlist[baseNodesIndex+23] = xyzmax[2];

	for(int i=0; i<8; i++){
		in.pointmtrlist[numInnerNodes+i] = size;		
		in.pointmarkerlist[numInnerNodes+i] = 1;
	}

	static int faceIndices[6][4]={
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{1, 2, 6, 5},
		{2, 3, 7, 6},
		{0, 3, 7, 4},
		{1, 0, 4, 5}};
	for(int i=0; i<zHandleSurface.triangles.size(); i++){
		in.facetmarkerlist[i] = 2;
        tetgenio::facet *f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = nullptr;
        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = zHandleSurface.triangles[i]->nodes[0]->index;
        p->vertexlist[1] = zHandleSurface.triangles[i]->nodes[1]->index;
        p->vertexlist[2] = zHandleSurface.triangles[i]->nodes[2]->index;	
	}

	//top
	in.facetmarkerlist[numInnerFacets] = 1;
	tetgenio::facet *f = &in.facetlist[numInnerFacets];
	f->numberofpolygons = 2;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 1;
	f->holelist = new double[3];
	tetgenio::polygon *p1 = &f->polygonlist[0];
	p1->numberofvertices = 4;
    p1->vertexlist = new int[p1->numberofvertices];	
	for(int i=0; i<4; i++){
		p1->vertexlist[i] = faceIndices[1][i]+numInnerNodes;
	}
	tetgenio::polygon *p2 = &f->polygonlist[1];
	p2->numberofvertices = topNodesList.size();
    p2->vertexlist = new int[p2->numberofvertices];	
	for(int i=0; i<p2->numberofvertices; i++){
		p2->vertexlist[i] = topNodesList[i]->index;
	}
	Vector3D holeCenter((oxyzmax[0]+oxyzmin[0])/2.0, (oxyzmax[1]+oxyzmin[1])/2.0, oxyzmax[2]);
	f->holelist[0] = holeCenter[0];
	f->holelist[1] = holeCenter[1];
	f->holelist[2] = holeCenter[2];

	//bottom	
	in.facetmarkerlist[numInnerFacets+1] = 1;
	f = &in.facetlist[numInnerFacets+1];
	f->numberofpolygons = 2;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 1;
	f->holelist = new double[3];
	p1 = &f->polygonlist[0];
	p1->numberofvertices = 4;
    p1->vertexlist = new int[p1->numberofvertices];	
	for(int i=0; i<4; i++){
		p1->vertexlist[i] = faceIndices[0][i]+numInnerNodes;
	}
	p2 = &f->polygonlist[1];
	p2->numberofvertices = bottomNodesList.size();
    p2->vertexlist = new int[p2->numberofvertices];	
	for(int i=0; i<p2->numberofvertices; i++){
		p2->vertexlist[i] = bottomNodesList[i]->index;
	}
	Vector3D holeCenter2((oxyzmax[0]+oxyzmin[0])/2.0, (oxyzmax[1]+oxyzmin[1])/2.0, oxyzmin[2]);
	f->holelist[0] = holeCenter2[0];
	f->holelist[1] = holeCenter2[1];
	f->holelist[2] = holeCenter2[2];	

	//left & right & back & front
	for(int i=2; i<6; i++){
		in.facetmarkerlist[numInnerFacets+i] = 1;
		tetgenio::facet *f = &in.facetlist[numInnerFacets+i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
        f->holelist = nullptr;
		tetgenio::polygon *p = &f->polygonlist[0];
		p->numberofvertices = 4;
		p->vertexlist = new int[p1->numberofvertices];	
		for(int j=0; j<4; j++){
			p->vertexlist[j] = faceIndices[i][j]+numInnerNodes;
		}	
	}

	Vector3D handleHole = Vector3D((xyzmax[0]+xyzmin[0])/2.0, (xyzmax[1]+xyzmin[1])/2.0, (xyzmin[2]+xyzmax[2])/2.0);
	in.holelist[0] = handleHole[0];
	in.holelist[1] = handleHole[1];
	in.holelist[2] = handleHole[2];

	char cmd[] = "pDqmfQ";
	tetrahedralize(cmd, &in, &out);
	meshOut.loadTETGENIO(out, true);
	int numTriFacets = out.numberoftrifaces;
	for(auto &n: meshOut.nodes){
		n->label = 3;
	}
	for (int i=0; i<numTriFacets; i++){
		if (out.trifacemarkerlist[i]==1){
			int base = i*3;
			meshOut.nodes[out.trifacelist[base]]->label = 1;
			meshOut.nodes[out.trifacelist[base+1]]->label = 1;
			meshOut.nodes[out.trifacelist[base+2]]->label = 1;
		}
	}
	for(auto &t: meshOut.tetrahedrons){
		t->label = 1;
	}
}

void parseZHandle(SurfaceMesh &zHandleSurface, double top, double bottom, 
	std::vector<std::array<double, 2>> &topEdgeNodes,	
    std::vector<std::array<double, 2>> &bottomEdgeNodes,
    std::vector<std::array<int, 2>> &topEdges, 
    std::vector<std::array<int, 2>> &bottomEdges){
	
	zHandleSurface.rebuildTriangleAdjacency();
	std::unordered_map<int, int> topNodeMap;
	std::unordered_map<int, int> bottomNodeMap;
	int topIndex = 0;
	int bottomIndex = 0;
	auto getTopNode =
	[&zHandleSurface, &topEdgeNodes, &topNodeMap, &topIndex]
	(int index){
		int rst = 0;
		if(topNodeMap.count(index)){
			rst = topNodeMap[index];
		}
		else{
			rst = topIndex++;
			topNodeMap[index]  = rst;
			topEdgeNodes.push_back({zHandleSurface.nodes[index]->pos[0], zHandleSurface.nodes[index]->pos[1]});
		}
		return rst;
	};
	auto getBottomNode =
	[&zHandleSurface, &bottomEdgeNodes, &bottomNodeMap, &bottomIndex]
	(int index){
		int rst = 0;
		if(bottomNodeMap.count(index)){
			rst = bottomNodeMap[index];
		}
		else{
			rst = bottomIndex++;
			bottomNodeMap[index]  = rst;
			bottomEdgeNodes.push_back({zHandleSurface.nodes[index]->pos[0], zHandleSurface.nodes[index]->pos[1]});
		}
		return rst;
	};
	for (auto &n: zHandleSurface.nodes){
		if (abs(n->pos[2]-top)<1e-13){
			n->edit = 1;
		}
		else if(abs(n->pos[2]-bottom)<1e-13){
			n->edit = 2;
		}
		else{
			n->edit = 0;
		}
	}
	std::vector<Triangle *> delTriangles;
	for (auto &t: zHandleSurface.triangles){
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
	for(auto &t: zHandleSurface.triangles){
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
	zHandleSurface.deleteTriangles(delTriangles);

}


void generateZHandleMeshV3(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet){
	Mesh goalMesh;
	Vector3D xyzmin;
	Vector3D xyzmax;
	Vector3D oxyzmin;
	Vector3D oxyzmax;
	tetgenio tetIn;
	std::vector<int> indexOf1;
	loadNodesWithLabel(tetIn, fileIn, xyzmax, xyzmin, oxyzmax, oxyzmin, indexOf1);

	if(indexOf1.size()!=8){
		std::cout << "Missing bounding box nodes!" << std::endl;
		exit(1);
	}

	Vector3D deltaPos(oxyzmax[0] - oxyzmin[0], oxyzmax[1] - oxyzmin[1], oxyzmax[2] - oxyzmin[2]);
	resetPoints(tetIn, oxyzmax+deltaPos, oxyzmin-deltaPos, indexOf1);
	Mesh innerMesh;
	SurfaceMesh interSurface;
	generateConvaxHullFromPointsIn3D(tetIn, oxyzmax, oxyzmin, innerMesh, interSurface, true);

	std::cout << innerMesh.nodes.size() << std::endl;
	for(auto &t: innerMesh.tetrahedrons){
		t->label = 0;
	}
	innerMesh.rebuildTetrahedronsAdjacency();
	// for(auto n: innerMesh.nodes){
	// 	n->label = 0;
	// }
	for(auto t: innerMesh.tetrahedrons){
		bool isContinue = false;
		for (int i=0; i<4; i++){
			if (t->nodes[i]->label != 0) {
				isContinue = true;
				break;
			}
		}
		if(isContinue){
			t->label = 1;
		}
	}
	
	for(auto t: innerMesh.tetrahedrons){
		if(t->label==1){
			for(auto v: t->nodes){
				if(v->label==0) v->label = 3;
			}
			for(int i=0; i<4; i++){
				if (t->adjacentTetrahedrons[i] && t->adjacentTetrahedrons[i]->label==0){
					SubTriangle f =  t->getSubTriangle(i);
					for(auto v: f.forms){
						v->label = 2;
					}
				}
				else if(!t->adjacentTetrahedrons[i]){
					SubTriangle f =  t->getSubTriangle(i);
					for(auto v: f.forms){
						if ((v->pos[2]>=(oxyzmax[2]-1e-6) || v->pos[2]<=(oxyzmin[2]+1e-6)) && (v->label==0 || v->label == 3)) v->label = 1;
					}
				}
			}
		}
		else if(t->label==0){
			for(int i=0; i<4; i++){
				if(!t->adjacentTetrahedrons[i]){
					SubTriangle f =  t->getSubTriangle(i);
					for(auto v: f.forms){
						v->label = 2;
					}
				}
			}
		}
	}

	innerMesh.exportVTK(fileOut+"_test.vtk");
	interSurface.exportVTK(fileOut+"_surf.vtk");
	Mesh out;
	parseZHandleV3(interSurface, xyzmax, xyzmin, oxyzmax, oxyzmin, size, out);
	innerMesh.mergeMesh(out, 1e-10);
	innerMesh.exportMESH(fileOut);
	innerMesh.exportVTK(fileOut+".vtk");

}








void generateZHandleMeshV2(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet){
	Mesh goalMesh;
	Vector3D xyzmin;
	Vector3D xyzmax;
	Vector3D oxyzmin;
	Vector3D oxyzmax;
	tetgenio tetIn;
	std::vector<int> indexOf1;
	loadNodesWithLabel(tetIn, fileIn, xyzmax, xyzmin, oxyzmax, oxyzmin, indexOf1);

	if(indexOf1.size()!=8){
		std::cout << "Missing bounding box nodes!" << std::endl;
		exit(1);
	}

	Vector3D deltaPos(oxyzmax[0] - oxyzmin[0], oxyzmax[1] - oxyzmin[1], oxyzmax[2] - oxyzmin[2]);
	resetPoints(tetIn, oxyzmax+deltaPos, oxyzmin-deltaPos, indexOf1);
	Mesh innerMesh;
	SurfaceMesh interSurface;
	generateConvaxHullFromPointsIn3D(tetIn, oxyzmax, oxyzmin, innerMesh, interSurface);
	std::cout << innerMesh.nodes.size() << std::endl;
	for(auto &t: innerMesh.tetrahedrons){
		t->label = 0;
	}
	innerMesh.rebuildTetrahedronsAdjacency();
	for(auto n: innerMesh.nodes){
		n->label = 0;
	}
	for(auto t: innerMesh.tetrahedrons){
		for (int i=0; i<4; i++){
			if (!t->adjacentTetrahedrons[i]){
				SubTriangle f = t->getSubTriangle(i);
				for(auto n: f.forms){
					n->label = 2;
				}
			}
		}
	}
	
	Mesh out;
	parseZHandleV3(interSurface, xyzmax, xyzmin, oxyzmax, oxyzmin, size, out);
	innerMesh.mergeMesh(out, 1e-10);
	innerMesh.exportMESH(fileOut);
	innerMesh.exportVTK(fileOut+".vtk");
}

void generateZHandleMesh(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet){
	Mesh goalMesh;
	Vector3D xyzmin;
	Vector3D xyzmax;
	Vector3D oxyzmin;
	Vector3D oxyzmax;
	tetgenio tetIn;
	std::vector<int> indexOf1;
	loadNodesWithLabel(tetIn, fileIn, xyzmax, xyzmin, oxyzmax, oxyzmin, indexOf1);

	if(indexOf1.size()!=8){
		std::cout << "Missing bounding box nodes!" << std::endl;
		exit(1);
	}

	Vector3D deltaPos(oxyzmax[0] - oxyzmin[0], oxyzmax[1] - oxyzmin[1], oxyzmax[2] - oxyzmin[2]);
	resetPoints(tetIn, oxyzmax+deltaPos, oxyzmin-deltaPos, indexOf1);
	Mesh innerMesh;
	SurfaceMesh interSurface;
	generateConvaxHullFromPointsIn3D(tetIn, oxyzmax, oxyzmin, innerMesh, interSurface);

	std::vector<std::array<int, 2>> bottomEdges;
	std::vector<std::array<double, 2>> bottomEdgeNodes;
	std::vector<std::array<int, 2>> topEdges;
	std::vector<std::array<double, 2>> topEdgeNodes;
	parseZHandle(interSurface, oxyzmax[2], oxyzmin[2], topEdgeNodes, bottomEdgeNodes, topEdges, bottomEdges);

	//Generate outer surface mesh	
	SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
	std::array<double, 2> oxymax({oxyzmax[0], oxyzmax[1]});
	std::array<double, 2> oxymin({oxyzmin[0], oxyzmin[1]});


	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, bottomEdgeNodes, bottomEdges);
	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, topEdgeNodes, topEdges);
	std::vector<std::array<double,2>> holes;
	holes.push_back({0.5*(oxymax[0]+oxymin[0]), 0.5*(oxymax[1]+oxymin[1])});
	triangulateio triBottom;
	triangulateio triTop;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triBottom);
	generateMeshInPlaneWithEdges(topEdgeNodes, topEdges, holes, size*size/2, triTop);

	faceBottom.projectTRIANGULATEIO(triBottom, PROJECTION_TYPE::XY_PLANE, xyzmin[2]);
	faceTop.projectTRIANGULATEIO(triTop, PROJECTION_TYPE::XY_PLANE, xyzmax[2]);
	deleteTRIANGULATEIOAllocatedArrays(triBottom);
	deleteTRIANGULATEIOAllocatedArrays(triTop);


	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[1], xyzmax[2]}, {xyzmin[1],xyzmin[2]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triFrontBack;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triFrontBack);
	faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmax[0]);
	faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmin[0]);
	deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

	bottomEdgeNodes.clear();
	bottomEdges.clear();
	holes.clear();
	generateRectangle({xyzmax[2], xyzmax[0]}, {xyzmin[2],xyzmin[0]}, size, bottomEdgeNodes, bottomEdges);
	triangulateio triLeftRight;	
	generateMeshInPlaneWithEdges(bottomEdgeNodes, bottomEdges, holes, size*size/2, triLeftRight);
	faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmax[1]);
	faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmin[1]);
	deleteTRIANGULATEIOAllocatedArrays(triLeftRight);





	interSurface.mergeSurfaceMesh(faceRight,1e-10);
	interSurface.mergeSurfaceMesh(faceLeft,1e-10);
	interSurface.mergeSurfaceMesh(faceFront,1e-10);
	interSurface.mergeSurfaceMesh(faceBack,1e-10);
	interSurface.mergeSurfaceMesh(faceTop,1e-10);	
	interSurface.mergeSurfaceMesh(faceBottom,1e-10);

	tetgenio tetLoop;
	interSurface.exportTETGENIO(tetLoop);
	tetLoop.numberofholes = 1;
	tetLoop.holelist = new double[3];

	interSurface.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/periodic/dislocation_out.vtk");

	return;
	// goalMesh.exportMESH(fileOut);
}


void generatePeriodicBoundaryConditionMesh(const std::string &fileIn, const std::string &fileOut, double size, bool beQuiet){
	Mesh goalMesh;
	Vector3D xyzmin;
	Vector3D xyzmax;
	Vector3D oxyzmin;
	Vector3D oxyzmax;
	tetgenio tetIn;
	std::vector<int> indexOf1;
	loadNodesWithLabel(tetIn, fileIn, xyzmax, xyzmin, oxyzmax, oxyzmin, indexOf1);

	if(indexOf1.size()!=8){
		std::cout << "Missing bounding box nodes!" << std::endl;
		exit(1);
	}

	Vector3D deltaPos(oxyzmax[0] - oxyzmin[0], oxyzmax[1] - oxyzmin[1], oxyzmax[2] - oxyzmin[2]);
	resetPoints(tetIn, oxyzmax+deltaPos, oxyzmin-deltaPos, indexOf1);
	Mesh innerMesh;
	SurfaceMesh interSurface;
	generateConvaxHullFromPointsIn3D(tetIn, oxyzmax, oxyzmin, innerMesh, interSurface);
	return;
	std::vector<std::array<double, 2>> planeNodes;
	for(int i=0; i<tetIn.numberofpoints; i++){
		if (abs(tetIn.pointlist[i*3+2] - xyzmin[2]) < 1e-13){
			std::array<double, 2> pos({tetIn.pointlist[i*3], tetIn.pointlist[i*3+1]});
			planeNodes.push_back(pos);
		}
	}


	std::array<double, 2> oxymax({oxyzmax[0], oxyzmax[1]});
	std::array<double, 2> oxymin({oxyzmin[0], oxyzmin[1]});
	std::vector<std::array<int, 2>> edges;
	std::vector<std::array<double, 2>> edgeNodes;
	generateConvaxHullFromPointsInPlane(planeNodes, oxymax, oxymin, edgeNodes, edges);

	generateRectangle({xyzmax[0],xyzmax[1]},{xyzmin[0], xyzmin[1]}, size, edgeNodes, edges);
	std::vector<std::array<double,2>> holes;
	holes.push_back({0.5*(oxymax[0]+oxymin[0]), 0.5*(oxymax[1]+oxymin[1])});
	triangulateio triOut;	
	generateMeshInPlaneWithEdges(edgeNodes, edges, holes, size*size/2, triOut);
	SurfaceMesh faceBottom, faceTop, faceFront, faceBack, faceLeft, faceRight;
	faceBottom.projectTRIANGULATEIO(triOut, PROJECTION_TYPE::XY_PLANE, xyzmin[2]);
	faceTop.projectTRIANGULATEIO(triOut, PROJECTION_TYPE::XY_PLANE, xyzmax[2]);
	deleteTRIANGULATEIOAllocatedArrays(triOut);


	edgeNodes.clear();
	edges.clear();
	holes.clear();
	generateRectangle({xyzmax[1], xyzmax[2]}, {xyzmin[1],xyzmin[2]}, size, edgeNodes, edges);
	triangulateio triFrontBack;	
	generateMeshInPlaneWithEdges(edgeNodes, edges, holes, size*size/2, triFrontBack);
	faceFront.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmax[0]);
	faceBack.projectTRIANGULATEIO(triFrontBack, PROJECTION_TYPE::YZ_PLANE, xyzmin[0]);
	deleteTRIANGULATEIOAllocatedArrays(triFrontBack);

	edgeNodes.clear();
	edges.clear();
	holes.clear();
	generateRectangle({xyzmax[2], xyzmax[0]}, {xyzmin[2],xyzmin[0]}, size, edgeNodes, edges);
	triangulateio triLeftRight;	
	generateMeshInPlaneWithEdges(edgeNodes, edges, holes, size*size/2, triLeftRight);
	faceRight.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmax[1]);
	faceLeft.projectTRIANGULATEIO(triLeftRight, PROJECTION_TYPE::ZX_PLANE, xyzmin[1]);
	deleteTRIANGULATEIOAllocatedArrays(triLeftRight);

	// goalMesh.exportMESH(fileOut);
}

void generateConvaxHullFromPointsIn3DRemoveHoles(tetgenio &tet, Mesh &goalMesh, SurfaceMesh &goalSurface){
	tetgenio tetout;
	char cmd[] = "Q";
	tetrahedralize(cmd, &tet, &tetout);

	goalMesh.loadTETGENIO(tetout);
	goalMesh.rebuildTetrahedronsAdjacency();


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
	for(auto &tet: goalMesh.tetrahedrons){
		tet->edit = 0;
		tetmesh.circumsphere( tet->nodes[0]->pos.data(),tet->nodes[1]->pos.data(),tet->nodes[2]->pos.data(),tet->nodes[3]->pos.data(), tmpCenter, &radius);
		tetRadius.push_back(radius);
	}
	std::vector<double> tmpTetRadius = tetRadius;
	std::sort(tmpTetRadius.begin(), tmpTetRadius.end());
	double eps = tmpTetRadius[tetRadius.size()*0.5];

	int removes=0;
	do{
		removes=0;
		for(int i = 0; i <goalMesh.tetrahedrons.size(); i++){
			if (goalMesh.tetrahedrons[i]->edit) continue;
			bool isBorder = false;
			for (int j = 0; j < 4; j++){
				if (goalMesh.tetrahedrons[i]->adjacentTetrahedrons[j]==nullptr || goalMesh.tetrahedrons[i]->adjacentTetrahedrons[j]->edit){
					isBorder = true;
					break;
				}
			}
			if (isBorder){
				if (tetRadius[i]>1.5*eps){
					goalMesh.tetrahedrons[i]->edit = 1;
					removes++;
				}
			}
		}
	}while(removes);

	for(auto &e: goalMesh.tetrahedrons){
		if(e->edit) continue;
		for(int i=0; i<4; i++){
			if (e->adjacentTetrahedrons[i]==nullptr || e->adjacentTetrahedrons[i]->edit){
				goalSurface.addTriangle(getNode(e->nodes[(i+1)%4]), getNode(e->nodes[(i+2)%4]), getNode(e->nodes[(i+3)%4]));
			}
		}
	}

	for(auto &n: oldNewNodes){
		goalSurface.nodes.push_back(n.second);
	}
	

	// for(int i = 0; i <goalMesh.tetrahedrons.size(); i++){
	// 	if (goalMesh.tetrahedrons[i]->edit) continue;
	// 	if (tetRadius[i]>1.5*eps){
	// 		goalMesh.tetrahedrons[i]->edit = 1;
	// 		removes++;
	// 	}
	// }

	for(int i=0; i<goalMesh.tetrahedrons.size();i++){
		if (goalMesh.tetrahedrons[i]->edit){
			delete goalMesh.tetrahedrons[i];
			goalMesh.tetrahedrons[i]=goalMesh.tetrahedrons.back();
			goalMesh.tetrahedrons.pop_back();
			//goalMesh.tetrahedrons.erase(goalMesh.tetrahedrons.begin()+i);
			i--;
		}
	}

	goalMesh.rebuildIndices();	


	// goalMesh.exportMESH("/home/kjfu/research/Mesher3DForSJTU/examples/two_holes2/inner3d.mesh");
	goalSurface.rebuildIndices();
}


void generateConvaxHullFromPointsIn3D(tetgenio &tet, Mesh &goalMesh, SurfaceMesh &goalSurface){
	tetgenio tetout;
	char cmd[] = "Q";
	tetrahedralize(cmd, &tet, &tetout);

	goalMesh.loadTETGENIO(tetout);
	goalMesh.rebuildTetrahedronsAdjacency();


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


	for(auto &e: goalMesh.tetrahedrons){

		for(int i=0; i<4; i++){

			if (e->adjacentTetrahedrons[i]==nullptr){
				Triangle *tri = new Triangle(getNode(e->nodes[(i+1)%4]), getNode(e->nodes[(i+2)%4]), getNode(e->nodes[(i+3)%4]));
				goalSurface.triangles.push_back(tri);
			}
		}
	}




	for(auto &n: oldNewNodes){
		goalSurface.nodes.push_back(n.second);
	}

	goalSurface.rebuildIndices();
}
void generateConvaxHullFromPointsIn3D(tetgenio &tet, Vector3D &oxyzmax, Vector3D &oxyzmin, Mesh &goalMesh, SurfaceMesh &goalSurface, bool withPointLabel){
	tetgenio tetout;
	char cmd[] = "Q";
	tetrahedralize(cmd, &tet, &tetout);
	// std::cout << tet.numberofpoints << "  " << tetout.numberofpoints << std::endl;
	goalMesh.loadTETGENIO(tetout, withPointLabel);
	goalMesh.rebuildTetrahedronsAdjacency();

	auto box3dContain=
    [&oxyzmax, &oxyzmin]
    (Vector3D &pos, double eps){
		bool rst=true;
        for(int i=0; i<3; i++){
            if(pos[i]>(oxyzmax[i]+eps) || (pos[i]<oxyzmin[i]-eps)){
                rst = false;
            }
        }
        return rst;
    };



	for(auto &n: goalMesh.nodes){
		n->edit = 0;
		if(!box3dContain(n->pos, 1e-8)){
			n->edit=1;
		}
	}

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
	for(auto &tet: goalMesh.tetrahedrons){
		tet->edit = 0;
		tetmesh.circumsphere( tet->nodes[0]->pos.data(),tet->nodes[1]->pos.data(),tet->nodes[2]->pos.data(),tet->nodes[3]->pos.data(), tmpCenter, &radius);
		tetRadius.push_back(radius);
		if (tet->nodes[0]->edit || tet->nodes[1]->edit || tet->nodes[2]->edit || tet->nodes[3]->edit) tet->edit = 1;
	}
	std::vector<double> tmpTetRadius = tetRadius;
	std::sort(tmpTetRadius.begin(), tmpTetRadius.end());
	double eps = tmpTetRadius[tetRadius.size()*0.5];

	int removes=0;
	do{
		removes=0;
		for(int i = 0; i <goalMesh.tetrahedrons.size(); i++){
			if (goalMesh.tetrahedrons[i]->edit) continue;
			bool isBorder = false;
			for (int j = 0; j < 4; j++){
				if (goalMesh.tetrahedrons[i]->adjacentTetrahedrons[j]==nullptr || goalMesh.tetrahedrons[i]->adjacentTetrahedrons[j]->edit){
					isBorder = true;
					break;
				}
			}
			if (isBorder){
				if (tetRadius[i]>2*eps){
					goalMesh.tetrahedrons[i]->edit = 1;
					removes++;
				}
			}
		}
	}while(removes);

	for(auto &e: goalMesh.tetrahedrons){
		if(e->edit) continue;
		for(int i=0; i<4; i++){
			if (e->adjacentTetrahedrons[i]==nullptr || e->adjacentTetrahedrons[i]->edit){
				Triangle *tri = new Triangle(getNode(e->nodes[(i+1)%4]), getNode(e->nodes[(i+2)%4]), getNode(e->nodes[(i+3)%4]));
				goalSurface.triangles.push_back(tri);
			}
		}
	}

	for(auto &n: oldNewNodes){
		goalSurface.nodes.push_back(n.second);
	}

	// for(auto &e: goalMesh.tetrahedrons){
	// 	e->edit=0;
	// 	for(int i=0; i<4; i++){
	// 		if(e->nodes[i]->edit==1){
	// 			e->edit=1;
	// 		}
	// 		if (e->nodes[i]->edit ==1 && e->nodes[(i+1)%4]->edit!=1 && e->nodes[(i+2)%4]->edit!=1 && e->nodes[(i+3)%4]->edit!=1 ){

	// 			TriangleElement *tri = new TriangleElement(getNode(e->nodes[(i+1)%4]), getNode(e->nodes[(i+2)%4]), getNode(e->nodes[(i+3)%4]));
	// 			goalSurface.triangles.push_back(tri);
	// 		}
	// 	}
	// }

	// for (int i=0; i<goalMesh.tetrahedrons.size(); i++){
	// 	if(goalMesh.tetrahedrons[i]->edit==1){
	// 		delete goalMesh.tetrahedrons[i];
	// 		goalMesh.tetrahedrons.erase(goalMesh.tetrahedrons.begin()+i);

	// 		i--;
	// 	}
	// }

	// for (int i=0; i<goalMesh.nodes.size(); i++){
	// 	if(goalMesh.nodes[i]->edit==1){
	// 		delete goalMesh.nodes[i];
	// 		goalMesh.nodes.erase(goalMesh.nodes.begin()+i);
	// 		i--;
	// 	}
	// }
	// goalMesh.rebuildIndices();
	for(int i=0; i<goalMesh.tetrahedrons.size();i++){
		if (goalMesh.tetrahedrons[i]->edit){
			delete goalMesh.tetrahedrons[i];
			goalMesh.tetrahedrons[i]=goalMesh.tetrahedrons.back();
			goalMesh.tetrahedrons.pop_back();
			//goalMesh.tetrahedrons.erase(goalMesh.tetrahedrons.begin()+i);
			i--;
		}
	}

	goalMesh.rebuildIndices();	

	// for(auto &n: oldNewNodes){
	// 	goalSurface.nodes.push_back(n.second);
	// }

	goalSurface.rebuildIndices();







	// goalMesh.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/PeriodicBoundaryConditionMesh/test.vtk");
	// goalSurface.exportVTK("/home/kjfu/research/Mesher3DForSJTU/examples/PeriodicBoundaryConditionMesh/test_s.vtk");
}

void generateConvaxHullFromPointsInPlane(std::vector<std::array<double, 2>> &planeNodes,
std::array<double, 2> &oxymax, std::array<double, 2> &oxymin,
std::vector<std::array<double,2>> &finalNodes, std::vector<std::array<int, 2>> &finalEdges){

	struct triangulateio triIn;
	struct triangulateio triOut;
	setNullToTRIANGULATEIO(triIn);
	setNullToTRIANGULATEIO(triOut);
	triIn.numberofpoints=planeNodes.size();
	triIn.numberofpointattributes =0;
	triIn.pointlist = new double[2*triIn.numberofpoints];
	for(int i=0; i<planeNodes.size(); i++){
		triIn.pointlist[i*2] = planeNodes[i][0];
		triIn.pointlist[i*2+1] = planeNodes[i][1];
	}
	char cmd[]="Q";
	triangulate(cmd, &triIn, &triOut, nullptr);
	SurfaceMesh aSurface;
	aSurface.projectTRIANGULATEIO(triOut, PROJECTION_TYPE::XY_PLANE, 0);
    
	deleteTRIANGULATEIOAllocatedArrays(triIn);
	deleteTRIANGULATEIOAllocatedArrays(triOut);

	auto box2dContain=
    [&oxymax, &oxymin]
    (Vector3D &pos, double eps){
		bool rst=true;
        for(int i=0; i<2; i++){
            if(pos[i]>(oxymax[i]+eps) || (pos[i]<oxymin[i]-eps)){
                rst = false;
            }
        }
        return rst;
    };

	for(auto &n: aSurface.nodes){
		n->edit=0;
		if (!box2dContain(n->pos, 1e-8)){
			n->edit=1;
		}
	}

	std::vector<std::array<Node*, 2>> tmpEdges;

	for(auto &e:aSurface.triangles){
		e->edit=0;
		for(int i=0; i<3; i++){
			if (e->nodes[i]->edit == 1 ){
				int next = (i+1)%3;
				int nextnext = (i+2)%3;
				if(e->nodes[next]->edit==0 && e->nodes[nextnext]->edit==0){
					tmpEdges.push_back({e->nodes[next], e->nodes[nextnext]});
				}
				break;
			}
		}
	}
	std::vector<Node*> outputNodes;
	for(auto &e:tmpEdges){

		e[0]->edit=2;
		e[1]->edit=2;
	}

	for(auto &n:aSurface.nodes){
		if(n->edit==2){
			n->index = outputNodes.size();
			outputNodes.push_back(n);
			finalNodes.push_back({n->pos[0], n->pos[1]});
		}
	}

	for(auto &e:tmpEdges){
		finalEdges.push_back({e[0]->index, e[1]->index});
	}

	// std::cout << finalEdges.size() << std::endl;
}



// void setNullToTRIANGULATEIO(triangulateio &io){
// 	io.edgelist = nullptr;
// 	io.edgemarkerlist = nullptr;
// 	io.holelist=nullptr;
// 	io.neighborlist=nullptr;
// 	io.normlist= nullptr;
// 	io.pointattributelist=nullptr;
// 	io.pointmarkerlist=nullptr;
// 	io.pointlist=nullptr;
// 	io.regionlist=nullptr;
// 	io.segmentlist=nullptr;
// 	io.segmentmarkerlist=nullptr;
// 	io.trianglearealist=nullptr;
// 	io.triangleattributelist=nullptr;
// 	io.trianglelist=nullptr;
// }

// void deleteTRIANGULATEIOAllocatedArrays(triangulateio &io){
// 	delete [] io.edgelist;
// 	delete [] io.edgemarkerlist;
// 	//delete [] io.holelist;
// 	delete [] io.neighborlist;
// 	delete [] io.normlist;
// 	delete [] io.pointattributelist;
// 	delete [] io.pointmarkerlist;
// 	delete [] io.pointlist;
// 	delete [] io.regionlist;
// 	delete [] io.segmentlist;
// 	delete [] io.segmentmarkerlist;
// 	delete [] io.trianglearealist;
// 	delete [] io.triangleattributelist;
// 	delete [] io.trianglelist;	
// }

void generateMeshInPlaneWithEdges(std::vector<std::array<double,2>> &planeNodes, std::vector<std::array<int, 2>> &edges, std::vector<std::array<double,2>> holes, double maxAreaSize,  triangulateio &triOut){
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


	std::string str = "pq30Qa"+std::to_string(maxAreaSize);
	char cmd[256];
	strcpy(cmd, str.c_str());
	triangulate(cmd, &triIn, &triOut, nullptr);
	deleteTRIANGULATEIOAllocatedArrays(triIn);
}


void generateRectangle(std::array<double, 2> maxPos, std::array<double,2> minPos, double size, std::vector<std::array<double,2>> &edgeNodes, std::vector<std::array<int, 2>> &edges){

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



void analysisQuality(const std::string &fileIn, const std::string &fileOut){
	Mesh mesh;
	mesh.loadMESH(fileIn);
	std::ofstream file(fileOut);
	file << "Quality" << std::endl;
	for(int i=0; i<mesh.tetrahedrons.size(); i++){
		Tetrahedron *tet =  mesh.tetrahedrons[i];
		file << calculateTetrahedronScaleQualityWith4Points_ISO(tet->nodes[0]->pos, tet->nodes[1]->pos, tet->nodes[2]->pos, tet->nodes[3]->pos) << std::endl; 
	}


}
}