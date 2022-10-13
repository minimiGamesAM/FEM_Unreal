// Fill out your copyright notice in the Description page of Project Settings.


#include "DynamicMesh.h"
#include "MeshDescriptionToDynamicMesh.h"
#include "StaticMeshAttributes.h"
#include "FemFunctions.h"
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
		
		FTriangleArray& triangles = MeshDescription->Triangles();
		std::vector<int> faces(triangles.Num() * 3, 0);
		
		for (const FTriangleID triId : triangles.GetElementIDs())
		{
			auto tri = MeshDescription->GetTriangleVertices(triId);
			
			faces[triId.GetValue() * 3]		= tri[0].GetValue();
			faces[triId.GetValue() * 3 + 1] = tri[1].GetValue();
			faces[triId.GetValue() * 3 + 2] = tri[2].GetValue();
		}

		////////////////////////////////////////////////////////////////////////

		int		numberOfPoints;
		double* pointlist;
		int		numberoftrifaces;
		int* trifacelist;
		int		numberoftetrahedra;
		int* tetrahedronlist;
		int* tet2facelist;

		UFemFunctions::runTetGen(points,
								 &faces[0],
								 faces.size(),
			                     VertexIDs.Num() * 3,
								 numberOfPoints,
								 pointlist,
								 numberoftrifaces,
								 trifacelist,
								 numberoftetrahedra,
								 tetrahedronlist,
								 tet2facelist,
								 "pqz-f-nn");

		MeshOut.Clear();


	}
}


// Called every frame
void ADynamicMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

