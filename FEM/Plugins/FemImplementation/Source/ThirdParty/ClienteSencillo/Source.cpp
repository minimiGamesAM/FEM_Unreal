#include<iostream>


#include "FemImpLibrary.h"

int main()
{
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

	std::cout << "Como unreal " << std::endl;
}