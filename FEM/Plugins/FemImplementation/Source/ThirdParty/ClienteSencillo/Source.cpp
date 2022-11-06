#include<iostream>

#include "FemImpLibrary.h"
#include <functional>
#include <algorithm>
#include <map>

//#include "TetGenInterface.h"

void testInvert(float* jac, const int dim)
{
		
	//float* jacCopy = new float[dim * dim];
	//
	//float* jacResult = new float[dim * dim];
	//
	//std::copy(jac, jac + dim * dim, jacCopy);
	//
	//float det = invert(jac, dim);
	//
	//matmul(jac, jacCopy, jacResult, dim, dim, dim);
	//	
	//for (int i = 0; i < dim * dim; ++i)
	//{
	//    std::cout << jacResult[i] << std::endl;
	//}
	//std::cout << "determinante " << det << std::endl;
	//
	//delete[] jacCopy;
	//delete[] jacResult;
}

void unitTest()
{
	//const int dim3 = 3;
	//const int dim4 = 4;
	//
	//float jac1[dim3 * dim3] = { 4, -2, 1,
	//							5, 0, 3,
	//							-1, 2, 6 };
	//
	//testInvert(jac1, dim3);
	//
	//std::cout << "*******************" << std::endl;
	//
	//float jac2[dim4 * dim4] = { 1, 1, 0, 3,
	//						   2, 1, -1, 1,
	//						   3, -1, -1, 2,
	//						   -1, 2, 3, -1 };
	//
	//testInvert(jac2, dim4);
	//
	//std::cout << "*******************" << std::endl;
	//
	//float jac3[dim3 * dim3] = { 0, 0, 1,
	//							1, 0, 0,
	//							1, 1, 0 };
	//
	//testInvert(jac3, dim3);
	//
	//std::cout << "*******************" << std::endl;
	//
	//float jac4[dim3 * dim3] = { 1, 3, -2,
	//							2, -3, 4,
	//							-1, 2, 5 };
	//
	//testInvert(jac4, dim3);
	//
	//std::cout << "*******************" << std::endl;
}

template<class T>
void prueba()
{
	const int dim = 3;
	//nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
	const int nodof = 3;
	//nn = total number of nodes in the problem
	const int nn = 8;


	const int verticesSize = nodof * nn;
	T vertices[verticesSize] = {

		T(0.0), T(1.0),	T( 0.0), T( 1.0), T(0.0), T(1.0), T( 0.0), T(1.0),
		T(0.0), T(0.0),	T( 0.0), T( 0.0), T(1.0), T(1.0), T( 1.0), T(1.0),
		T(0.0), T(0.0), T(-1.0), T(-1.0), T(0.0), T(0.0), T(-1.0), T(-1.0) };

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
	
	int nf[nodof * nn] = { 0, 1, 0, 1, 0, 1, 0, 1,
						   0, 0, 0, 0, 1, 1, 1, 1,
						   1, 1, 0, 0, 1, 1, 0, 0 };

	// number of nodes per element
	int nod = tetsSize / nbTets;

	int id = FEM_Factory<T>::create(dim, nodof, nbTets, nod, 1, "tetrahedron");
	FEM_Factory<T>::init(id, vertices, tets, nf, nn);

	std::vector<T> verticesBuffer(nn * 3, T(0.0));

	for (int i = 1; i <= 20; ++i)
	{
		FEM_Factory<T>::update(id, T(1.0), &verticesBuffer[0]);
	}

	//node 5
	//std::cout << verticesBuffer[12] << " " << verticesBuffer[13] << " " << verticesBuffer[14] << std::endl;
}

template<class T>
void prueba2D()
{
	const int dim = 2;
	//nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
	const int nodof = 2;
	//nn = total number of nodes in the problem
	const int nn = 18;

	// nels = nb elements
	const int nels = 3;

	const int nip = 9;

	
	T vertices[] = {
		T(0.0000E+00),
		T(0.0000E+00),
		T(0.0000E+00),
		T(0.6667E+00),
		T(0.6667E+00),
		T(0.1333E+01),
		T(0.1333E+01),
		T(0.1333E+01),
		T(0.2000E+01),
		T(0.2000E+01),
		T(0.2667E+01),
		T(0.2667E+01),
		T(0.2667E+01),
		T(0.3333E+01),
		T(0.3333E+01),
		T(0.4000E+01),
		T(0.4000E+01),
		T(0.4000E+01),
		T(0.0000E+00),
		T(-0.5000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.5000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.5000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.1000E+01),
		T(0.0000E+00),
		T(-0.5000E+00),
		T(-0.1000E+01)
	};

	int g_num[] = {
		 3,
		 2,
		 1,
		 4,
		 6,
		 7,
		 8,
		 5,
		 8,
		 7,
		 6,
		 9,
		11,
		12,
		13,
		10,
		13,
		12,
		11,
		14,
		16,
		17,
		18,
		15
	};

	std::for_each(std::begin(g_num), std::end(g_num), [](int& i) { i = i - 1; });

	// number of nodes per element
	int nod = std::distance(std::begin(g_num), std::end(g_num)) / nels;
	int id = FEM_Factory<T>::create(dim, nodof, nels, nod, nip, "quadrilateral");

	int nf[] = {
		0,
		0,
		0,
		0,
		0,
		0,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1,
		1
	};
	
	std::vector<int> tempnf(nodof * nn, 0);
	
	for (int k = 0; k < nn; ++k)
	{
		for (int l = 0; l < nodof; ++l)
		{
			tempnf[k + l * nn] = nf[l + k * nodof];
		}
	}

	FEM_Factory<T>::init(id, vertices, g_num, &tempnf[0], nn);

	///
	std::vector<int> nodesLoaded;
	nodesLoaded.push_back(18);

	std::vector<T> val;
	val.push_back(T(0.0));
	val.push_back(T(1.0));

	FEM_Factory<T>::loadedNodes(id, &nodesLoaded[0], nodesLoaded.size(), &val[0]);

	/////
	for (int i = 1; i <= 20; ++i)
	{
		FEM_Factory<T>::update(id, T(1.0), nullptr);
	}
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
		

	//prueba<double>();

	prueba2D<double>();

}