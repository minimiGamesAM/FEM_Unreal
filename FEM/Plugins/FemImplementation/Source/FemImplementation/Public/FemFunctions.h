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
	static void runFem(float* verticesBuffer, int bufferSize);
};
