// Fill out your copyright notice in the Description page of Project Settings.


#include "FemFunctions.h"

THIRD_PARTY_INCLUDES_START
#include "FemImpLibrary/FemImpLibrary.h"
THIRD_PARTY_INCLUDES_END

void UFemFunctions::runFem(float* verticesBuffer, int bufferSize)
{

	float determinante = basicTest(verticesBuffer, bufferSize);
	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Determinante %f"), determinante));
}