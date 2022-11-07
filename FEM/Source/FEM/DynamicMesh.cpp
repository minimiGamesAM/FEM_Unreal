// Fill out your copyright notice in the Description page of Project Settings.


#include "DynamicMesh.h"
#include "MeshDescriptionToDynamicMesh.h"
#include "DynamicMeshToMeshDescription.h"
#include "StaticMeshAttributes.h"
#include "FemFunctions.h"
#include <vector>

#include "Generators/GridBoxMeshGenerator.h"

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

	void UpdateStaticMeshFromDynamicMesh(
		UStaticMesh* StaticMesh,
		const FDynamicMesh3* Mesh)
	{
		FMeshDescription MeshDescription;
		FStaticMeshAttributes StaticMeshAttributes(MeshDescription);
		StaticMeshAttributes.Register();

		FDynamicMeshToMeshDescription Converter;
		Converter.Convert(Mesh, MeshDescription);

		// todo: vertex color support

		//UStaticMesh* StaticMesh = NewObject<UStaticMesh>(Component);
		//FName MaterialSlotName = StaticMesh->AddMaterial(MyMaterial);

		// Build the static mesh render data, one FMeshDescription* per LOD.
		TArray<const FMeshDescription*> MeshDescriptionPtrs;
		MeshDescriptionPtrs.Emplace(&MeshDescription);
		StaticMesh->BuildFromMeshDescriptions(MeshDescriptionPtrs);
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

	TArray<UActorComponent*> comps;

	this->GetComponents(comps);
	
	TArray<FColor>	VertexColors;

	for (int i = 0; i < comps.Num(); ++i) //Because there may be more components
	{
		UProceduralMeshComponent* thisComp = Cast<UProceduralMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisComp)
		{
			auto meshSection = thisComp->GetProcMeshSection(0);
			
			auto& data = meshSection->ProcVertexBuffer;
				
			mVerticesBuffer.resize(data.Num() * 3);
						
			for (int j = 0; j < data.Num(); ++j)
			{
				//position
				FVector& pos = data[j].Position;
								
				mVerticesBuffer[j * 3] = pos[0];
				mVerticesBuffer[j * 3 + 1] = pos[1];
				mVerticesBuffer[j * 3 + 2] = pos[2];

				//color
				FColor& color = data[j].Color;
				VertexColors.Add(color);
			}
		}
	}

	// init FEM
	const int dim = 3;
	//nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
	const int nodof = 3;
	//nn = total number of nodes in the problem
	int nn = mVerticesBuffer.size() / 3;

	std::vector<int> nf(nodof * nn, 1);

	for (int i = 0; i < VertexColors.Num(); ++i)
	{
		if (VertexColors[i] == FColor(255, 0, 0, 255))
		{
			nf[i] = 0;
			nf[i + nn] = 0;
			nf[i + 2 * nn] = 0;
		}
	}

	std::vector<float> tempVerticesBuffer(nodof * nn, 0.0f);

	for (int k = 0; k < nn; ++k)
	{
		for (int l = 0; l < nodof; ++l)
		{
			tempVerticesBuffer[k + l * nn] = mVerticesBuffer[l + k * nodof];
		}
	}

	// number of nodes per element
	int nod = 4;

	mIdAlgoFEM = UFemFunctions::create(dim, nodof, mTetsBuffer.size() / 4, nod, 1, "tetrahedron");
	UFemFunctions::init(mIdAlgoFEM, &tempVerticesBuffer[0], &mTetsBuffer[0], &nf[0], nn);

	///////////
	std::vector<int> nodesLoaded;
	nodesLoaded.push_back(4);

	std::vector<float> val;
	val.push_back(float(0.33));
	val.push_back(float(0.33));
	val.push_back(float(0.33));

	UFemFunctions::loadedNodes(mIdAlgoFEM, &nodesLoaded[0], nodesLoaded.size(), &val[0]);
	///////////

	for (int i = 0; i < 20; ++i)
	{
		UFemFunctions::update(mIdAlgoFEM, 1.0, &mVerticesBuffer[0]);
	}

}

void ADynamicMesh::PostInitializeComponents()
{
	Super::PostInitializeComponents();
	OnMeshGenerationSettingsModified();
}

////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////

void ADynamicMesh::OnMeshGenerationSettingsModified()
{
	//EditMesh([this](FDynamicMesh3& MeshToUpdate) {
		RegenerateSourceMesh(SourceMesh);
	//	});
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
		//auto m = &StaticMesh->GetRenderData()->LODResources[0];
		//auto sec = m->Sections[0];

		TMap<FVector, FColor> VertexColorData;
		StaticMesh->GetVertexColorData(VertexColorData);
		
		//FColorVertexBuffer* colorBuffer = &StaticMesh->GetRenderData()->LODResources[0].VertexBuffers.ColorVertexBuffer;
		
		//colorBuffer->Init(colorBuffer->GetNumVertices(), true);
		//TArray<FColor> colors; 
		//colorBuffer->GetVertexColors(colors);
		
		FMeshDescription* MeshDescription = StaticMesh->GetMeshDescription(0);
		
		const FVertexArray& VertexIDs = MeshDescription->Vertices();
		
		//auto nbVertices = colorBuffer->GetNumVertices();

		//GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Pos %i, %i"), nbVertices, VertexIDs.Num()));

		//if (nbVertices == VertexIDs.Num())
		//{

		//}

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
		
		UFemFunctions::runTetGenio(	points,
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
		
		ProceduralMesh = Cast<UProceduralMeshComponent>(GetDefaultSubobjectByName(TEXT("ProcMesh")));

		if (ProceduralMesh)
		{
			ProceduralMesh->ClearAllMeshSections();

			//ProceduralMesh->CreateMeshSection(0, )

			auto defaultPos = FVector(0.0f, 0.0f, 0.0f);
			auto defaultColor = FColor(0, 0, 0, 0);

			TArray<FVector>				Vertices(&defaultPos, numberOfPoints);
			TArray<int32>				Triangles;
			TArray<FVector>				Normals;
			TArray<FVector2D>			UVs;
			TArray<FColor>		        VertexColors(&defaultColor, numberOfPoints);
			TArray<FProcMeshTangent>	Tangents;
			
			//vertex
			for (int i = 0; i < numberOfPoints; i++)
			{
				//MeshOut.AppendVertex(FVector(pointlist[i * 3], pointlist[i * 3] + 1, pointlist[i * 3] + 2));

				FVector p(	pointlist[i * 3 + 0],
							pointlist[i * 3 + 1],
							pointlist[i * 3 + 2]);

				UVs.Add(FVector2D(0.0, 0.0));
				UVs.Add(FVector2D(0.0, 1.0));
				UVs.Add(FVector2D(1.0, 0.0));

				Vertices[i] = p;
			}
			
			//colors
			for (int i = 0; i < Vertices.Num(); i++)
			{
				FVector& p = Vertices[i];

				FColor* itt = VertexColorData.Find(p);

				if (itt != nullptr)
				{
					VertexColors[i] = *itt;
				}
			}
			

			//triangles
			for (int tetIdx = 0; tetIdx < numberoftetrahedra; tetIdx++)
			{
				for (int i = 0; i < 4; ++i)
				{
					int idx = tet2facelist[(tetIdx * 4) + i];

					Triangles.Add(trifacelist[idx * 3 + 0]);
					Triangles.Add(trifacelist[idx * 3 + 1]);
					Triangles.Add(trifacelist[idx * 3 + 2]);

					//auto tri = FIndex3i(trifacelist[idx * 3 + 0],
					//	trifacelist[idx * 3 + 1],
					//	trifacelist[idx * 3 + 2]);

					//MeshOut.AppendTriangle(tri);

					//MeshOut.SetVertexUV(tri[0], FVector2f(0.0, 0.0));
					//MeshOut.SetVertexUV(tri[1], FVector2f(0.0, 1.0));
					//MeshOut.SetVertexUV(tri[2], FVector2f(1.0, 0.0));
				}
			}

			///////////////////////////
			mTetsBuffer.clear();
			mTetsBuffer.resize(numberoftetrahedra * 4);

			for (int tetIdx = 0; tetIdx < numberoftetrahedra; tetIdx++)
			{
				mTetsBuffer[4 * tetIdx] = tetrahedronlist[4 * tetIdx];
				mTetsBuffer[4 * tetIdx + 1] = tetrahedronlist[4 * tetIdx + 2];
				mTetsBuffer[4 * tetIdx + 2] = tetrahedronlist[4 * tetIdx + 1];
				mTetsBuffer[4 * tetIdx + 3] = tetrahedronlist[4 * tetIdx + 3];
			}
			///////////////////////////

			ProceduralMesh->CreateMeshSection(0, Vertices, Triangles, Normals, UVs, VertexColors, Tangents, false);
		
			UMaterialInterface* UseMaterial = UMaterial::GetDefaultMaterial(MD_Surface);
			
			if (Material != nullptr)
			{
				UseMaterial = Material;
			}

			ProceduralMesh->SetMaterial(0, UseMaterial);
		}
	}
}

// Called every frame
void ADynamicMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
			
	FVector defaultVec(0.0f, 0.0f, 0.0f);
	TArray<FVector>	Vertices(&defaultVec, mVerticesBuffer.size() / 3);
		
	// update array with FEM

	UFemFunctions::update(mIdAlgoFEM, 1.0, &mVerticesBuffer[0]);

	//for (int i = 0; i < mVerticesBuffer.size(); ++i)
	//{
	//	mVerticesBuffer[i] += 0.2 * t;
	//}
	//

	for (int i = 0; i < Vertices.Num(); ++i)
	{
		Vertices[i][0] = mVerticesBuffer[i * 3];
		Vertices[i][1] = mVerticesBuffer[i * 3 + 1];
		Vertices[i][2] = mVerticesBuffer[i * 3 + 2];
	}

	TArray<int32>				Triangles;
	TArray<FVector>				Normals;
	TArray<FVector2D>			UVs;
	TArray<FColor>		        VertexColors;
	TArray<FProcMeshTangent>	Tangents;

	ProceduralMesh->UpdateMeshSection(0, Vertices, Normals, UVs, VertexColors, Tangents);

}

