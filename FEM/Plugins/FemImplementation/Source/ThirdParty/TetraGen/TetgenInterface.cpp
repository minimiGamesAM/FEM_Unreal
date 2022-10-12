#include "TetGenInterface.h"

#include "tetgen/tetgen.h"
#include <string>

namespace tetSpace
{
	struct TetraedrosStruct
	{
		int   m_numberOfPoints;
		REAL *m_pointlist;
		int   m_numberoftrifaces;
		int  *m_trifacelist;
		int   m_numberoftetrahedra;
		int  *m_tetrahedronlist;
		int  *m_tet2facelist;

		TetraedrosStruct()
			: m_numberOfPoints(0),
			  m_pointlist(nullptr),
			  m_numberoftrifaces(0),
			  m_trifacelist(nullptr),
			  m_numberoftetrahedra(0),
			  m_tetrahedronlist(nullptr),
			  m_tet2facelist(nullptr)
		{
			//m_pointlist = new REAL[m_numberOfPoints];
			//trifacelist = new int[m_numberoftrifaces];
			//tetrahedronlist = new int[m_numberoftetrahedra];
		}
	};

	TetraedrosStruct tetraedrosInstancia;
}

void runTetGen2( float* points,
				 int* faces,
				 int facesSize,
				 int pointsSizes,
				 char* switches )
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

	tetSpace::tetraedrosInstancia.m_numberOfPoints = out.numberofpoints;
	tetSpace::tetraedrosInstancia.m_numberoftrifaces = out.numberoftrifaces;
	tetSpace::tetraedrosInstancia.m_numberoftetrahedra = out.numberoftetrahedra;

	////////////////////////////

	delete[] tetSpace::tetraedrosInstancia.m_pointlist;
	delete[] tetSpace::tetraedrosInstancia.m_trifacelist;
	delete[] tetSpace::tetraedrosInstancia.m_tetrahedronlist;
	delete[] tetSpace::tetraedrosInstancia.m_tet2facelist;

	tetSpace::tetraedrosInstancia.m_pointlist = new REAL[3 * out.numberofpoints];
	tetSpace::tetraedrosInstancia.m_trifacelist = new int[3 * out.numberoftrifaces];
	tetSpace::tetraedrosInstancia.m_tetrahedronlist = new int[4 * out.numberoftetrahedra];
	tetSpace::tetraedrosInstancia.m_tet2facelist = new int[4 * out.numberoftetrahedra];

	std::memcpy(tetSpace::tetraedrosInstancia.m_pointlist,
		out.pointlist, 3 * out.numberofpoints * sizeof(REAL));

	std::memcpy(tetSpace::tetraedrosInstancia.m_trifacelist,
		out.trifacelist, 3 * out.numberoftrifaces * sizeof(int));

	std::memcpy(tetSpace::tetraedrosInstancia.m_tetrahedronlist,
		out.tetrahedronlist, 4 * out.numberoftetrahedra * sizeof(int));

	std::memcpy(tetSpace::tetraedrosInstancia.m_tet2facelist,
		out.tet2facelist, 4 * out.numberoftetrahedra * sizeof(int));

}

int getNumberOfPoints()
{
	return tetSpace::tetraedrosInstancia.m_numberOfPoints;
}

double getPoint(int idx)
{
	return tetSpace::tetraedrosInstancia.m_pointlist[idx];
}

int getNumberOfTrifaces()
{
	return tetSpace::tetraedrosInstancia.m_numberoftrifaces;
}

int getTrifacet(int idx)
{
	return tetSpace::tetraedrosInstancia.m_trifacelist[idx];
}

int getNumberOfTets()
{
	return tetSpace::tetraedrosInstancia.m_numberoftetrahedra;
}

int getTet(int idx)
{
	return tetSpace::tetraedrosInstancia.m_tetrahedronlist[idx];
}

int getTet2facelist(int tetIdx, int idx)
{
	return tetSpace::tetraedrosInstancia.m_tet2facelist[tetIdx + idx];
}