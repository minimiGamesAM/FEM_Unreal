// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralMeshComponent.h"
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

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	//std::vector<float>			mVerticesBuffer;
};
