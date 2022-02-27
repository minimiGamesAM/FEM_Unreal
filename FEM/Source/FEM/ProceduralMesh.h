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
};
