#include<iostream>

#include "FemImpLibrary.h"
#include <functional>
#include <algorithm>
#include <map>

#include "TetGenInterface.h"

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

	int id = FEM_Factory<T>::create(dim, nodof, nbTets);
	FEM_Factory<T>::init(id, vertices, tets, nf, nn);

	std::vector<T> verticesBuffer(nn * 3, T(0.0));

	for (int i = 1; i <= 20; ++i)
	{
		FEM_Factory<T>::update(id, T(1.0), &verticesBuffer[0]);
	}

	//node 5
	//std::cout << verticesBuffer[12] << " " << verticesBuffer[13] << " " << verticesBuffer[14] << std::endl;
}

void testWithTetMesh()
{
	const int dim = 3;
	const int nodof = 3;

	char file[] = "verificacion.ply";
	char switches[] = "pqz-f-nn";

	runTetGen(file, switches);

	int nn = getNumberOfPoints();

	std::vector<float> mVerticesBuffer(nn * nodof, 0.0f);

	std::cout << "Vertices" << std::endl;

	for (int i = 0; i < nn; i++)
	{
		mVerticesBuffer[i * 3] = getPoint(i * 3);
		mVerticesBuffer[i * 3 + 1] = getPoint(i * 3 + 1);
		mVerticesBuffer[i * 3 + 2] = getPoint(i * 3 + 2);

		std::cout << i << " " << mVerticesBuffer[i * 3] << " " << mVerticesBuffer[i * 3 + 1] << " " << mVerticesBuffer[i * 3 + 2] << std::endl;
	}
		
	std::vector<int> mTetsBuffer(getNumberOfTets() * 4, 0);

	std::cout << "Tets" << std::endl;

	for (int tetIdx = 0; tetIdx < getNumberOfTets(); tetIdx++)
	{
		mTetsBuffer[4 * tetIdx] = getTet(4 * tetIdx);
		mTetsBuffer[4 * tetIdx + 1] = getTet(4 * tetIdx + 2);
		mTetsBuffer[4 * tetIdx + 2] = getTet(4 * tetIdx + 1);
		mTetsBuffer[4 * tetIdx + 3] = getTet(4 * tetIdx + 3);

		std::cout << tetIdx << " " << mTetsBuffer[4 * tetIdx] << " "
									   << mTetsBuffer[4 * tetIdx + 1] << " " 
										<< mTetsBuffer[4 * tetIdx + 2] << " "
										<< mTetsBuffer[4 * tetIdx + 3]  << std::endl;
	}

	//mTetsBuffer[0] = 5;
	//mTetsBuffer[1] = 0;
	//mTetsBuffer[2] = 1;
	//mTetsBuffer[3] = 7;
	//////////

	std::vector<int> nf( { 0, 1, 0, 1, 0, 1, 0, 1,
						   0, 0, 0, 0, 1, 1, 1, 1,
						   1, 1, 0, 0, 1, 1, 0, 0 });


	/////////////////transpose////////////////////////////////////
	std::vector<float> tempVertices(nodof * nn, 0);

	for (int k = 0; k < nn; ++k)
	{
		for (int l = 0; l < nodof; ++l)
		{
			tempVertices[k + l * nn] = mVerticesBuffer[l + k * nodof];
		}
	}


	int mIdAlgoFEM = FEM_Factory<float>::create(dim, nodof, getNumberOfTets());
	FEM_Factory<float>::init(mIdAlgoFEM, &tempVertices[0], &mTetsBuffer[0], &nf[0], nn);

	for (int i = 1; i <= 200; ++i)
	{
		FEM_Factory<float>::update(mIdAlgoFEM, float(0.25), &mVerticesBuffer[0]);
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

	testWithTetMesh();

	//prueba<double>();

}