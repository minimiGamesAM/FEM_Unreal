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

void runTetGen(char* file, char* switches)
{
	tetgenio in, out;
	tetgenio addin, bgmin;

	in.load_ply(const_cast<char*>(std::string(file).c_str()));
			
	// Output the PLC to files 'barin.node' and 'barin.poly'.
	//in.save_nodes("barin");
	//in.save_poly("barin");

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).

	//char switches[] = "pq1.414a0.1";

	tetrahedralize(switches, &in, &out, &addin, &bgmin);

	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	//out.save_nodes("barout");
	//out.save_elements("barout");
	//out.save_faces("barout");

	tetSpace::tetraedrosInstancia.m_numberOfPoints		= out.numberofpoints;
	tetSpace::tetraedrosInstancia.m_numberoftrifaces	= out.numberoftrifaces;
	tetSpace::tetraedrosInstancia.m_numberoftetrahedra	= out.numberoftetrahedra;

	////////////////////////////

	delete[] tetSpace::tetraedrosInstancia.m_pointlist;
	delete[] tetSpace::tetraedrosInstancia.m_trifacelist;
	delete[] tetSpace::tetraedrosInstancia.m_tetrahedronlist;
	delete[] tetSpace::tetraedrosInstancia.m_tet2facelist;

	tetSpace::tetraedrosInstancia.m_pointlist		= new REAL[3 * out.numberofpoints];
	tetSpace::tetraedrosInstancia.m_trifacelist		= new int[3 * out.numberoftrifaces];
	tetSpace::tetraedrosInstancia.m_tetrahedronlist = new int[4 * out.numberoftetrahedra];
	tetSpace::tetraedrosInstancia.m_tet2facelist	= new int[4 * out.numberoftetrahedra];
	
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