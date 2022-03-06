#include "TetGenInterface.h"

#include "tetgen/tetgen.h"
#include <string>

tetgenio in, out;
tetgenio addin, bgmin;

void runTetGen(char* file, char* switches)
{
	in.clean_memory();
	in.initialize();

	out.clean_memory();
	out.initialize();

	in.load_stl(const_cast<char*>(std::string(file).c_str()));
			
	// Output the PLC to files 'barin.node' and 'barin.poly'.
	in.save_nodes("barin");
	in.save_poly("barin");

	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).

	//char switches[] = "pq1.414a0.1";

	tetrahedralize(switches, &in, &out, &addin, &bgmin);

	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	out.save_nodes("barout");
	out.save_elements("barout");
	out.save_faces("barout");
}

int getNumberOfPoints()
{
	return out.numberofpoints;
}

double getPoint(int idx)
{
	return out.pointlist[idx];
}

int getNumberOfTrifaces()
{
	return out.numberoftrifaces;
}

int getTrifacet(int idx)
{
	return out.trifacelist[idx];
}

int getNumberOfTets()
{
	return out.numberoftetrahedra;
}

int getTet(int idx)
{
	return out.tetrahedronlist[idx];
}

int getTet2facelist(int tetIdx, int idx)
{
	return out.tet2facelist[tetIdx][idx];
}
