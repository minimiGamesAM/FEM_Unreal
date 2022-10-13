// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "FemFunctions.generated.h"

/**
 * 
 */
UCLASS()
class FEMIMPLEMENTATION_API UFemFunctions : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
	
private:

public:
	static int create(int ndim, int nodof, int nels);
	static void init(int id, float* g_coord, int* g_num, int* in_nf, int in_nn);
	static void update(int id, float dt, float* verticesBuffer);

	static void runTetGen(const float* points,
		const int*	faces,
		const int	facesSize,
		const int	pointsSizes,
		int&		numberOfPoints,
		double*&	pointlist,
		int&		numberoftrifaces,
		int*&		trifacelist,
		int&		numberoftetrahedra,
		int*&		tetrahedronlist,
		int*&		tet2facelist,
		char*		switches);


	static void runFem(float* verticesBuffer, int verticesBufferSize, int* tetsBuffer, int tetsBufferSize);
};
