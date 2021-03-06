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

	const int verticesSize = 18;
	float vertices[verticesSize] = {
		0, 0, 1,
		0, 1, 1,
		0, 0, 0,
		0, 1, 0,
		1, 0, 0,
		1, 1, 0
	};
	
	const int tetsSize = 3 * 4;
	int tets[tetsSize] = {
		0, 4, 5, 2,
		1, 0, 5, 2,
		1, 2, 5, 3
	};
	
	int tetId = 0;
	elemStiffnessMatrix(vertices, &tets[tetId * 4]);
	
	const int nodof = 3;
	const int nbNodes = verticesSize / 3;
	int* nf = new int[nodof * verticesSize / 3];

	int count = 1;
	std::map<int, bool> nodeInactif;

	nodeInactif[2] = true;
	nodeInactif[3] = true;
	nodeInactif[4] = true;
	
	for (int i = 0; i < nbNodes; ++i)
	{
		for (int j = 0; j < nodof; ++j)
		{
			if (nodeInactif.count(i))
			{
				nf[i * nodof + j] = 0;
			}
			else
			{
				nf[i * nodof + j] = count;
				count++;
			}
		}
	}

	for (int i = 0; i < nbNodes * nodof; ++i)
	{
		std::cout << "nf " << nf[i] << std::endl;
	}
}