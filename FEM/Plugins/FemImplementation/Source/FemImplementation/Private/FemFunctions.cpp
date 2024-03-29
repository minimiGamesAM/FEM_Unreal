// Fill out your copyright notice in the Description page of Project Settings.


#include "FemFunctions.h"

THIRD_PARTY_INCLUDES_START
#include "FemImpLibrary/FemImpLibrary.h"
#include "TetraGen/TetGenInterface.h"
THIRD_PARTY_INCLUDES_END

int UFemFunctions::create(int ndim, int nodof, int nels, int nod, int nip, const char* element)
{
	return FEM_Factory<float>::create(ndim, nodof, nels, nod, nip, element);
}

void UFemFunctions::loadedNodes(int id, int* nodes, int loaded_nodes, float* vals)
{
	FEM_Factory<float>::loadedNodes(id, nodes, loaded_nodes, vals);
}

void UFemFunctions::setGravityAcceleration(int id, float gravAcc)
{
	FEM_Factory<float>::setGravityAcceleration(id, gravAcc);
}

void UFemFunctions::setGravityDirection(int id, float* gravDir)
{
	FEM_Factory<float>::setGravityDirection(id, gravDir);
}

void UFemFunctions::setDamping(int id, const float fk, const float fm)
{
	FEM_Factory<float>::setDamping(id, fk, fm);
}

void UFemFunctions::setMaterialParams(int id, const float e, const float v, const float gamma)
{
	FEM_Factory<float>::setMaterialParams(id, e, v, gamma);
}

void UFemFunctions::init(int id, float* g_coord, int* g_num, int* in_nf, int in_nn)
{
	FEM_Factory<float>::init(id, g_coord, g_num, in_nf, in_nn);
}

long long UFemFunctions::update(int id, float dt, float* verticesBuffer)
{
	return FEM_Factory<float>::update(id, dt, verticesBuffer);
}

void UFemFunctions::runFem(float* verticesBuffer, int verticesBufferSize, int* tetsBuffer, int tetsBufferSize)
{
	float determinante = basicTest(verticesBuffer, verticesBufferSize, tetsBuffer, tetsBufferSize);
	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Determinante %f"), determinante));
}

void UFemFunctions::runTetGenio(const float* points,
	const int* faces,
	const int	facesSize,
	const int	pointsSizes,
	int&			numberOfPoints,
	double*& pointlist,
	int&			numberoftrifaces,
	int*& trifacelist,
	int&	numberoftetrahedra,
	int*& tetrahedronlist,
	int*& tet2facelist,
	char* switches)
{
	runTetGen(points,
		faces,
		facesSize,
		pointsSizes,
		numberOfPoints,
		pointlist,
		numberoftrifaces,
		trifacelist,
		numberoftetrahedra,
		tetrahedronlist,
		tet2facelist,
		switches);
}