// Fill out your copyright notice in the Description page of Project Settings.


#include "DynamicMesh.h"
#include "TetGenFunctionLibrary.h"
#include "MeshDescriptionToDynamicMesh.h"
#include "StaticMeshAttributes.h"
#include <vector>

namespace
{
	void UpdateDynamicMeshFromStaticMesh(
		const UStaticMesh* StaticMesh,
		FDynamicMesh3& Mesh)
	{
		FMeshDescription* MeshDescription = StaticMesh->GetMeshDescription(0);

		FMeshDescriptionToDynamicMesh converter;
		converter.Convert(MeshDescription, Mesh);
	}
}
// Sets default values
ADynamicMesh::ADynamicMesh()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

}

// Called when the game starts or when spawned
void ADynamicMesh::BeginPlay()
{
	Super::BeginPlay();
	
}

void ADynamicMesh::PostLoad()
{
	Super::PostLoad();
	OnMeshGenerationSettingsModified();
}

void ADynamicMesh::PostActorCreated()
{
	Super::PostActorCreated();
	OnMeshGenerationSettingsModified();
}

void ADynamicMesh::EditMesh(TFunctionRef<void(FDynamicMesh3&)> EditFunc)
{
	EditFunc(SourceMesh);

	// update spatial data structures
	//if (bEnableSpatialQueries || bEnableInsideQueries)
	//{
	//	MeshAABBTree.Build();
	//	if (bEnableInsideQueries)
	//	{
	//		FastWinding->Build();
	//	}
	//}

	OnMeshEditedInternal();

}
void ADynamicMesh::OnMeshGenerationSettingsModified()
{
	EditMesh([this](FDynamicMesh3& MeshToUpdate) {
		RegenerateSourceMesh(MeshToUpdate);
		});
}

void ADynamicMesh::OnMeshEditedInternal()
{
	if (MeshComponent == nullptr)
	{
		MeshComponent = Cast<UStaticMeshComponent>(GetDefaultSubobjectByName(TEXT("Mesh")));

		if (MeshComponent)
		{
			StaticMesh = MeshComponent->GetStaticMesh();
		}
	}

	OnMeshModified.Broadcast(this);
}

void ADynamicMesh::RegenerateSourceMesh(FDynamicMesh3& MeshOut)
{
	if (MeshComponent == nullptr)
	{
		MeshComponent = Cast<UStaticMeshComponent>(GetDefaultSubobjectByName(TEXT("Mesh")));

		if (MeshComponent)
		{
			StaticMesh = MeshComponent->GetStaticMesh();
		}
	}

	if (MeshComponent && StaticMesh)
	{

		FMeshDescription* MeshDescription = StaticMesh->GetMeshDescription(0);

		const FVertexArray& VertexIDs = MeshDescription->Vertices();

		TVertexAttributesConstRef<FVector> VertexPositions =
			MeshDescription->VertexAttributes().GetAttributesRef<FVector>(MeshAttribute::Vertex::Position);

		float* points = new float[VertexIDs.Num() * 3];

		for (const FVertexID vId : VertexIDs.GetElementIDs())
		{
			FVector p = VertexPositions.Get(vId);

			points[vId.GetValue() * 3] = p.X;
			points[vId.GetValue() * 3 + 1] = p.Y;
			points[vId.GetValue() * 3 + 2] = p.Z;
		}
		
		//const FPolygonArray& Polygons = MeshDescription->Polygons();
		//std::vector<std::vector<int>> facets;
		//
		//for (const FPolygonID PolygonID : Polygons.GetElementIDs())
		//{
		//	int32 PolygonGroupID = MeshDescription->GetPolygonPolygonGroup(PolygonID).GetValue();
		//	const TArray<FTriangleID>& TriangleIDs = MeshDescription->GetPolygonTriangleIDs(PolygonID);
		//	int NumTriangles = TriangleIDs.Num();
		//	facets.emplace_back();// [PolygonID.GetValue()] .resize(facets.size() + NumTriangles);
		//
		//	for (int TriIdx = 0; TriIdx < NumTriangles; ++TriIdx)
		//	{
		//		auto InstanceTri = MeshDescription->GetTriangleVertices(TriangleIDs[TriIdx]);
		//		facets.back().push_back(InstanceTri[0].GetValue());
		//		facets.back().push_back(InstanceTri[1].GetValue());
		//		facets.back().push_back(InstanceTri[2].GetValue());
		//
		//		GEngine->AddOnScreenDebugMessage(-1, 100.0f, FColor::Yellow, FString::Printf(TEXT("facet %i %i %i"), InstanceTri[0].GetValue(), InstanceTri[1].GetValue(), InstanceTri[2].GetValue()));
		//	}
		//
		//	GEngine->AddOnScreenDebugMessage(-1, 100.0f, FColor::Yellow, FString::Printf(TEXT("***************")));
		//}

		FTriangleArray& triangles = MeshDescription->Triangles();
		std::vector<int> faces(triangles.Num() * 3, 0);
		
		for (const FTriangleID triId : triangles.GetElementIDs())
		{
			auto tri = MeshDescription->GetTriangleVertices(triId);
			
			faces[triId.GetValue() * 3]		= tri[0].GetValue();
			faces[triId.GetValue() * 3 + 1] = tri[1].GetValue();
			faces[triId.GetValue() * 3 + 2] = tri[2].GetValue();
		}

		UTetGenFunctionLibrary::RunTetGen2(points, &faces[0], faces.size(), VertexIDs.Num() * 3, nullptr);
	}


}


// Called every frame
void ADynamicMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

