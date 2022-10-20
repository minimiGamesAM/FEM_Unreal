#include "TetGenInterface.h"

#include "tetgen/tetgen.h"
#include <string>


void runTetGen( const float* points,
				 const int* faces,
				 const int facesSize,
				 const int pointsSizes,
				 int&       numberOfPoints,
				 double*&  pointlist,
				 int&       numberoftrifaces,
				 int*&     trifacelist,
				 int&       numberoftetrahedra,
				 int*&     tetrahedronlist,
				 int*&     tet2facelist,
				 char*	   switches )
{
	tetgenio in, out;
	tetgenio::facet* f;
	tetgenio::polygon* p;
	
	in.firstnumber = 0;
	in.numberofpoints = pointsSizes / 3;
	
	in.pointlist = new REAL[in.numberofpoints * 3];

	for (int i = 0; i < in.numberofpoints; ++i)
	{
		in.pointlist[i * 3 + 0] = points[i * 3 + 0];
		in.pointlist[i * 3 + 1] = points[i * 3 + 1];
		in.pointlist[i * 3 + 2] = points[i * 3 + 2];
	}

	in.numberoffacets = facesSize / 3;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = nullptr;

	for (int i = 0; i < in.numberoffacets; ++i)
	{
		f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;

		for (int j = 0; j < f->numberofpolygons; ++j)
		{
			//int face_index = faces[];

			p = &f->polygonlist[j];
			p->numberofvertices = 3;
			p->vertexlist = new int[p->numberofvertices];

			p->vertexlist[0] = faces[i * 3 + 0];
			p->vertexlist[1] = faces[i * 3 + 1];
			p->vertexlist[2] = faces[i * 3 + 2];
		}
	}

	tetrahedralize(switches, &in, &out);

	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	//out.save_nodes("barout");
	//out.save_elements("barout");
	//out.save_faces("barout");

	numberOfPoints = out.numberofpoints;
	numberoftrifaces = out.numberoftrifaces;
	numberoftetrahedra = out.numberoftetrahedra;

	////////////////////////////
	pointlist		= new double[3 * out.numberofpoints];
	trifacelist	= new int[3 * out.numberoftrifaces];
	tetrahedronlist = new int[4 * out.numberoftetrahedra];
	tet2facelist	= new int[4 * out.numberoftetrahedra];

	std::memcpy(pointlist,
		out.pointlist, 3 * out.numberofpoints * sizeof(REAL));

	std::memcpy(trifacelist,
		out.trifacelist, 3 * out.numberoftrifaces * sizeof(int));

	std::memcpy(tetrahedronlist,
		out.tetrahedronlist, 4 * out.numberoftetrahedra * sizeof(int));

	std::memcpy(tet2facelist,
		out.tet2facelist, 4 * out.numberoftetrahedra * sizeof(int));

}