// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"

#include <vector>
#include "Deformable.generated.h"

UCLASS()
class FEM_API ADeformable : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ADeformable();

	UPROPERTY(EditAnywhere)
	UProceduralMeshComponent* ProceduralMeshComp = nullptr;

	UPROPERTY(EditAnywhere)
	UStaticMeshComponent* StaticMeshComp = nullptr;

	UPROPERTY(EditAnywhere, Category = MaterialOptions)
	UMaterialInterface* Material = nullptr;

	UPROPERTY(EditAnywhere, Category = FEM_MaterialParams)
		float e = 10000.0f;

	UPROPERTY(EditAnywhere, Category = FEM_MaterialParams)
		float v = 0.3f;

	UPROPERTY(EditAnywhere, Category = FEM_MaterialParams)
		float gamma = 1.0f;

	UPROPERTY(EditAnywhere, Category = FEM_MaterialParams)
		float fk = 0.01f;

	UPROPERTY(EditAnywhere, Category = FEM_MaterialParams)
		float fm = 0.0f;

	UPROPERTY(EditAnywhere, Category = FEM_WorldParams)
		float gravity = 980.0f;

	UPROPERTY(EditAnywhere, Category = FEM_WorldParams)
		FVector gravityDir;
protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	//std::vector<float>			mVerticesBuffer;
private:
	void buildTetMesh();

private:
	std::vector<float>			mVerticesBuffer;
	std::vector<int>			mTetsBuffer;

	int							mIdAlgoFEM = -1;
};
