// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "DynamicMesh3.h"
#include "Components/StaticMeshComponent.h"
#include "DynamicMesh.generated.h"

UCLASS()
class FEM_API ADynamicMesh : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ADynamicMesh();

	UPROPERTY(VisibleAnywhere)
	UStaticMeshComponent* MeshComponent = nullptr;

	UPROPERTY(VisibleAnywhere)
	UStaticMesh* StaticMesh = nullptr;

	/**
	 * This delegate is broadcast whenever the internal SourceMesh is updated
	 */
	DECLARE_MULTICAST_DELEGATE_OneParam(FOnMeshModified, ADynamicMesh*);
	FOnMeshModified OnMeshModified;

	virtual void EditMesh(TFunctionRef<void(FDynamicMesh3&)> EditFunc);


protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

	virtual void PostLoad() override;
	virtual void PostActorCreated() override;

	virtual void OnMeshEditedInternal();

	/** Called whenever the initial Source mesh needs to be regenerated / re-imported. Calls EditMesh() to do so. */
	virtual void OnMeshGenerationSettingsModified();

	/** Called to generate or import a new source mesh. Override this to provide your own generated mesh. */
	virtual void RegenerateSourceMesh(FDynamicMesh3& MeshOut);

	FDynamicMesh3 SourceMesh;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

};
