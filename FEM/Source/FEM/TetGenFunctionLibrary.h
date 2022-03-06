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

	UFUNCTION(BlueprintCallable)
	static int getNumberOfPoints();

	UFUNCTION(BlueprintCallable)
	static float getPoint(int idx);

	UFUNCTION(BlueprintCallable)
	static int getNumberOfTrifaces();

	UFUNCTION(BlueprintCallable)
	static int getTrifacet(int idx);

	UFUNCTION(BlueprintCallable)
	static int getNumberOfTets();

	UFUNCTION(BlueprintCallable)
	static int getTet(int idx);

	UFUNCTION(BlueprintCallable)
	static int getTet2facelist(int tetIdx, int idx);
};
