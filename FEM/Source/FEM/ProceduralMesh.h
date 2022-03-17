// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Containers/Array.h"
#include "ProceduralMeshComponent.h"
#include "ProceduralMesh.generated.h"


UCLASS()
class FEM_API AProceduralMesh : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AProceduralMesh();

	UFUNCTION(BlueprintCallable)
	void runTetragenio();

	UFUNCTION(BlueprintCallable)
	TArray<FVector>& getVerticess();

	UFUNCTION(BlueprintCallable)
	TArray<int32>& getTriangulos();

	UFUNCTION(BlueprintCallable)
	TArray<FVector2D>& getUVs();

	//UFUNCTION(BlueprintCallable)
	//TArray<FVector>& getVertices();
//

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

private:

	TArray<FVector>				mVertices;
	TArray<int32>				mTriangles;
	TArray<FVector>				mNormals;
	TArray<FVector2D>			mUVs;
	TArray<FColor>		        mVertexColors;
	TArray<FProcMeshTangent>	mTangents;

	FProcMeshSection*			mMeshSection;

	float t = 0.0f;
};
