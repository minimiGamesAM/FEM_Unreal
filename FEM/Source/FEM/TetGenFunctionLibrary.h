// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "TetGenFunctionLibrary.generated.h"

/**
 * 
 */
UCLASS()
class FEM_API UTetGenFunctionLibrary : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()

private:
	static void* DllHandle;
	static bool LoadDllHandle();

public:
	UFUNCTION(BlueprintCallable)
	static void RunTetGen();
};
