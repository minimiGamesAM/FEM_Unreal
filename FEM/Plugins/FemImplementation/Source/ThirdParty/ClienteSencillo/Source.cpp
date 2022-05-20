#include<iostream>


#include "FemImpLibrary.h"

int main()
{
	float* buffer = new float[2];

	buffer[0] = 2;
	buffer[1] = 3;

	float det = basicTest(buffer, 2);

	std::cout << "Como unreal " << det << std::endl;
}