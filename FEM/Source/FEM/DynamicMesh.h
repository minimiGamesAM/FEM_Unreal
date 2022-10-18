// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "DynamicMesh3.h"
#include "Components/StaticMeshComponent.h"
#include "ProceduralMeshComponent.h"
#include "DynamicMesh.generated.h"

UCLASS()
class FEM_API ADynamicMesh : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ADynamicMesh();

	UPROPERTY(EditAnywhere)
	UStaticMeshComponent* MeshComponent = nullptr;

	UPROPERTY(EditAnywhere)
	UProceduralMeshComponent* ProceduralMesh = nullptr;

	UPROPERTY(EditAnywhere, Category = MaterialOptions)
	UMaterialInterface* Material = nullptr;

	UPROPERTY(Transient)
	UStaticMesh* StaticMesh = nullptr;

	/**
	 * This delegate is broadcast whenever the internal SourceMesh is updated
	 */
	//DECLARE_MULTICAST_DELEGATE_OneParam(FOnMeshModified, ADynamicMesh*);
	//FOnMeshModified OnMeshModified;

	//virtual void EditMesh(TFunctionRef<void(FDynamicMesh3&)> EditFunc);
	
protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

	//////////////////////////////////////////////
	//mutually exclusive
	virtual void PostLoad() override;
	virtual void PostActorCreated() override;
	//////////////////////////////////////////////

	/** Allow actors to initialize themselves on the C++ side after all of their components have been initialized, only called during gameplay */
	virtual void PostInitializeComponents();
		
	/** Called whenever the initial Source mesh needs to be regenerated / re-imported. Calls EditMesh() to do so. */
	//virtual void OnMeshGenerationSettingsModified();

	/** Called to generate or import a new source mesh. Override this to provide your own generated mesh. */
	virtual void RegenerateSourceMesh(FDynamicMesh3& MeshOut);

	void OnMeshGenerationSettingsModified();

	FDynamicMesh3 SourceMesh;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

private:
	TArray<FVector>				mVertices;

	float t = 0.0f;

};
