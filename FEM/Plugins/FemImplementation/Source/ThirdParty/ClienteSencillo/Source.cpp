#include<iostream>

#include "FemImpLibrary.h"
#include <functional>
#include <algorithm>
#include <map>

void testInvert(float* jac, const int dim)
{
		
	float* jacCopy = new float[dim * dim];

	float* jacResult = new float[dim * dim];

	std::copy(jac, jac + dim * dim, jacCopy);

	float det = invert(jac, dim);

	matmul(jac, jacCopy, jacResult, dim, dim, dim);
		
	for (int i = 0; i < dim * dim; ++i)
	{
	    std::cout << jacResult[i] << std::endl;
	}
	std::cout << "determinante " << det << std::endl;

	delete[] jacCopy;
	delete[] jacResult;
	//float jac[4 * 4] = { 1, 1, 0, 3,
	//                     2, 1, -1, 1,
	//                     3, -1, -1, 2,
	//                     -1, 2, 3, -1 };
	//
}

void unitTest()
{
	const int dim3 = 3;
	const int dim4 = 4;

	float jac1[dim3 * dim3] = { 4, -2, 1,
								5, 0, 3,
								-1, 2, 6 };

	testInvert(jac1, dim3);

	std::cout << "*******************" << std::endl;

	float jac2[dim4 * dim4] = { 1, 1, 0, 3,
							   2, 1, -1, 1,
							   3, -1, -1, 2,
							   -1, 2, 3, -1 };

	testInvert(jac2, dim4);

	std::cout << "*******************" << std::endl;

	float jac3[dim3 * dim3] = { 0, 0, 1,
								1, 0, 0,
								1, 1, 0 };

	testInvert(jac3, dim3);

	std::cout << "*******************" << std::endl;

	float jac4[dim3 * dim3] = { 1, 3, -2,
								2, -3, 4,
								-1, 2, 5 };

	testInvert(jac4, dim3);

	std::cout << "*******************" << std::endl;
}

int main()
{
	FEM_Factory factory = FEM_Factory();
	
	//unitTest();
	
	//float* buffer = new float[2];
	//
	//buffer[0] = 2;
	//buffer[1] = 3;
	//
	//int* bufferTets = new int[2];
	//
	//bufferTets[0] = 2;
	//bufferTets[1] = 3;

	//float det = basicTest(buffer, 2, bufferTets, 2);

	const int verticesSize = 24;
	float vertices[verticesSize] = {

		0.0f, 1.0f,	 0.0f,  1.0f, 0.0f, 1.0f,  0.0f, 1.0f,
		0.0f, 0.0f,	 0.0f,  0.0f, 1.0f, 1.0f,  1.0f, 1.0f,
		0.0f, 0.0f, -1.0f, -1.0f, 0.0f, 0.0f, -1.0f, -1.0f };
		
	const int nbTets = 6;
	const int tetsSize = nbTets * 4;
	int tets[tetsSize] = {
		1, 3, 4, 7,
		1, 4, 2, 7,
		1, 2, 5, 7,
		6, 4, 8, 7,
		6, 2, 4, 7,
		6, 5, 2, 7
	};

	std::for_each(tets, tets + tetsSize, [](int& i) { i = i - 1; });
	
	//loads in kN (node 1, 2, 5, 6)
	// 1  0.0  0.0 - 0.1667   2  0.0  0.0 - 0.3333
	// 5  0.0  0.0 - 0.3333   6  0.0  0.0 - 0.1667

	const int dim = 3;
	const int loadsNodeSize = 4;

	float loads[loadsNodeSize * dim] = {	0.0f,  0.0f, -0.1667f,
											0.0f,  0.0f, -0.3333f,
											0.0f,  0.0f, -0.3333f,
											0.0f,  0.0f, -0.1667f
	};

	int loads_nodes_ids[loadsNodeSize] = { 1,
									 2,
									 5,
									 6
	};
		
	elemStiffnessMatrix(vertices, tets, loads, nbTets, loads_nodes_ids, loadsNodeSize);// &tets[tetId * 4]);
	
	//nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
	const int nodof = 3;
	//nn = total number of nodes in the problem
	const int nn = 8;
	
	int nf[nodof * nn] = { 0, 1, 0, 1, 0, 1, 0, 1,
							0, 0, 0, 0, 1, 1, 1, 1,
							1, 1, 0, 0, 1, 1, 0, 0 };

	
	factory.create(dim, nodof, nbTets);
	factory.init(vertices, tets, nf, nn);
	factory.update();

}