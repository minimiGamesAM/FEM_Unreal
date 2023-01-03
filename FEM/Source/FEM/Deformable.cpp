// Fill out your copyright notice in the Description page of Project Settings.


#include "Deformable.h"
#include "MeshDescriptionToDynamicMesh.h"
#include "StaticMeshAttributes.h"

// Sets default values
ADeformable::ADeformable()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

}

// Called when the game starts or when spawned
void ADeformable::BeginPlay()
{
	Super::BeginPlay();

	TArray<UActorComponent*> comps;

	this->GetComponents(comps);
	
	for (int i = 0; i < comps.Num(); ++i) //Because there may be more components
	{
		UProceduralMeshComponent* thisCompP = Cast<UProceduralMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisCompP)
		{
			ProceduralMeshComp = thisCompP;
		}

		UStaticMeshComponent* thisCompS = Cast<UStaticMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisCompS)
		{
			StaticMeshComp = thisCompS;
		}
	}

	if (ProceduralMeshComp && StaticMeshComp)
	{
		//static mesh
		TMap<FVector, FColor> VertexColorData;
		StaticMeshComp->GetStaticMesh()->GetVertexColorData(VertexColorData);

		FMeshDescription* MeshDescription = StaticMeshComp->GetStaticMesh()->GetMeshDescription(0);
		const FVertexArray& VertexIDs = MeshDescription->Vertices();

		TVertexAttributesConstRef<FVector> VertexPositions =
		MeshDescription->VertexAttributes().GetAttributesRef<FVector>(MeshAttribute::Vertex::Position);

		// procedural mesh
		ProceduralMeshComp->ClearAllMeshSections();

		auto defaultPos = FVector(0.0f, 0.0f, 0.0f);
		auto defaultColor = FColor(0, 0, 0, 255);

		int		numberOfPoints(VertexPositions.GetNumElements());

		TArray<FVector>				Vertices;// (&defaultPos, numberOfPoints);
		TArray<int32>				Triangles;
		TArray<FVector>				Normals;
		TArray<FVector2D>			UVs;
		TArray<FColor>		        VertexColors(&defaultColor, numberOfPoints);
		TArray<FProcMeshTangent>	Tangents;

		//////////////
		for (const FVertexID vId : VertexIDs.GetElementIDs())
		{
			FVector p = VertexPositions.Get(vId);
			Vertices.Push(p);
		}

		FTriangleArray& triangles = MeshDescription->Triangles();

		for (const FTriangleID triId : triangles.GetElementIDs())
		{
			auto tri = MeshDescription->GetTriangleVertices(triId);

			Triangles.Push(tri[0].GetValue());
			Triangles.Push(tri[1].GetValue());
			Triangles.Push(tri[2].GetValue());
		}
		
		//////////////

		ProceduralMeshComp->CreateMeshSection(0, Vertices, Triangles, Normals, UVs, VertexColors, Tangents, false);

		UMaterialInterface* UseMaterial = UMaterial::GetDefaultMaterial(MD_Surface);

		if (Material != nullptr)
		{
			UseMaterial = Material;
		}

		ProceduralMeshComp->SetMaterial(0, UseMaterial);

		//auto meshSection = ProceduralMeshComp->GetProcMeshSection(0);

		//auto& data = meshSection->ProcVertexBuffer;

	}


}

// Called every frame
void ADeformable::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

