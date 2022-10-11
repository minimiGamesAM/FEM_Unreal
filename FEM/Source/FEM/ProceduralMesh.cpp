// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralMesh.h"
#include "TetGenFunctionLibrary.h"
#include "FemFunctions.h"
#include <algorithm>
#include <vector>
#include <iterator>

//https://nerivec.github.io/old-ue4-wiki/pages/procedural-mesh-component-in-cgetting-started.html
// without the "Pluging part"
// 1) Add in MyProject.Build.cs "ProceduralMeshComponent" at the end of PublicDependencyModuleNames.AddRange
// 2) Add in AdditionalDependencies in .uproject file "ProceduralMeshComponent"
// 3) To fix errors with Visual Studio IntelliSense you need to right-click MyProject.uproject and re-generate Visual Studio project files. In Visual Studio 2017, open "Solution Explorer" and open the "Game" folder, right-click on the first line, which should be the root of your solution, select: "Rescan Solution".


// Sets default values
AProceduralMesh::AProceduralMesh()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	
}

void AProceduralMesh::runTetragenio()
{
	mVertices.Empty();
	mTriangles.Empty();
	mNormals.Empty();
	mUVs.Empty();
	mVertexColors.Empty();
	mTangents.Empty();
	mTetsIds.Empty();
		
	//UTetGenFunctionLibrary::RunTetGen();
	int nbPoints = UTetGenFunctionLibrary::getNumberOfPoints();

	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("NB POINTS %i"), nbPoints));

	for (int i = 0; i < nbPoints; i++)
	{
		FVector p(UTetGenFunctionLibrary::getPoint(i * 3),
			UTetGenFunctionLibrary::getPoint(i * 3 + 1),
			UTetGenFunctionLibrary::getPoint(i * 3 + 2));

		mUVs.Add(FVector2D(0.0, 0.0));
		mUVs.Add(FVector2D(0.0, 1.0));
		mUVs.Add(FVector2D(1.0, 0.0));

		//GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Pos %f, %f, %f"), p.X, p.Y, p.Z));
		//to be remove this scale, should do a proper blender model with the right scale 
		mVertices.Add(p);
	}

	for (int tetIdx = 0; tetIdx < UTetGenFunctionLibrary::getNumberOfTets(); tetIdx++)
	{
		for (int i = 0; i < 4; ++i)
		{
			int idx = UTetGenFunctionLibrary::getTet2facelist(tetIdx * 4, i);
	
			mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(idx * 3 + 0));
			mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(idx * 3 + 1));
			mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(idx * 3 + 2));
		}
	}

	//for (int i = 0; i < UTetGenFunctionLibrary::getNumberOfTrifaces(); ++i)
	//{
	//	mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(i * 3));
	//	mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(i * 3 + 1));
	//	mTriangles.Add(UTetGenFunctionLibrary::getTrifacet(i * 3 + 2));
	//}
	
	mTetsIds.SetNum(UTetGenFunctionLibrary::getNumberOfTets() * 4);

	for (int tetIdx = 0; tetIdx < UTetGenFunctionLibrary::getNumberOfTets(); tetIdx++)
	{
		mTetsIds[4 * tetIdx] = UTetGenFunctionLibrary::getTet(4 * tetIdx);
		mTetsIds[4 * tetIdx + 1] = UTetGenFunctionLibrary::getTet(4 * tetIdx + 2);
		mTetsIds[4 * tetIdx + 2] = UTetGenFunctionLibrary::getTet(4 * tetIdx + 1);
		mTetsIds[4 * tetIdx + 3] = UTetGenFunctionLibrary::getTet(4 * tetIdx + 3);
	}

	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("NB Tets %i "), UTetGenFunctionLibrary::getNumberOfTets()));
	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("NB Facets %i "), UTetGenFunctionLibrary::getNumberOfTrifaces()));

}

// Called when the game starts or when spawned
void AProceduralMesh::BeginPlay()
{
	Super::BeginPlay();

	//runTetragenio();

	//mVertices.Empty();
	//mTriangles.Empty();
	//mNormals.Empty();
	//mUVs.Empty();
	//mVertexColors.Empty();
	//mTangents.Empty();
		
	TArray<UActorComponent*> comps;

	this->GetComponents(comps);

	for (int i = 0; i < comps.Num(); ++i) //Because there may be more components
	{
		UProceduralMeshComponent* thisComp = Cast<UProceduralMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisComp)
		{
			mMeshComponent = thisComp;
			//auto meshSection = thisComp->GetProcMeshSection(0);
			//
			//auto& data = meshSection->ProcVertexBuffer; // .Position;
			//	
			//mVertices.SetNum(data.Num());
			//
			//if (!mVerticesBuffer)
			//{
			//	mVerticesBuffer = new float[data.Num() * 3];
			//}
			//
			//for (int j = 0; j < data.Num(); ++j)
			//{
			//	mVertices[j] = data[j].Position;
			//
			//	mVerticesBuffer[j * 3] = mVertices[j][0];
			//	mVerticesBuffer[j * 3 + 1] = mVertices[j][1];
			//	mVerticesBuffer[j * 3 + 2] = mVertices[j][2];
			//}
		}

		/////////////////////////////////
		if (!mTetsBuffer)
		{
			mTetsBuffer = new int[mTetsIds.Num()];
		}

		for (int j = 0; j < mTetsIds.Num(); ++j)
		{
			mTetsBuffer[j] = mTetsIds[j];
		}

		/////////////////////////////////
		if (mVerticesBuffer.empty())
		{
			mVerticesBuffer.resize(mVertices.Num() * 3, 0.0f);
		}
		
		for (int j = 0; j < mVertices.Num(); ++j)
		{
			mVerticesBuffer[j * 3] = mVertices[j][0];
			mVerticesBuffer[j * 3 + 1] = mVertices[j][1];
			mVerticesBuffer[j * 3 + 2] = mVertices[j][2];
		}

		if (mIdAlgoFEM < 0)
		{
			const int dim = 3;
			//nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
			const int nodof = 3;

			//nn = total number of nodes in the problem
			const int nn = 8;

			nf = new int[nodof * nn];

			int temp_nf[nodof * nn] = { 0, 1, 0, 1, 0, 1, 0, 1,
										0, 0, 0, 0, 1, 1, 1, 1,
										1, 1, 0, 0, 1, 1, 0, 0 };

			std::copy(temp_nf, temp_nf + nodof * nn, nf);

			/////////////////transpose////////////////////////////////////
			std::vector<float> tempVertices(nodof * mVertices.Num(), 0);

			for (int k = 0; k < mVertices.Num(); ++k)
			{
				for (int l = 0; l < nodof; ++l)
				{
					tempVertices[k + l * mVertices.Num()] = mVerticesBuffer[l + k * nodof];
				}
			}

			//std::copy(std::begin(tempVertices), std::end(tempVertices), mVerticesBuffer);
			//////////////////////////////////////////////////////////////
						
			//mIdAlgoFEM = UFemFunctions::create(dim, nodof, mTetsIds.Num() / 4);
			//UFemFunctions::init(mIdAlgoFEM, &tempVertices[0], &mTetsBuffer[0], nf, nn);
		}
	}
}

TArray<FVector>& AProceduralMesh::getVerticess()
{
	return mVertices;
}

TArray<int32>& AProceduralMesh::getTriangulos()
{
	return mTriangles;
}

TArray<FVector2D>& AProceduralMesh::getUVs()
{
	return mUVs;
}

//void AProceduralMesh::storeTetsIds()
//{
//	mTetsIds.Empty();
//
//	mTetsIds.SetNum(UTetGenFunctionLibrary::getNumberOfTets() * 4);
//
//	for (int tetIdx = 0; tetIdx < UTetGenFunctionLibrary::getNumberOfTets(); tetIdx++)
//	{
//		mTetsIds[tetIdx] = UTetGenFunctionLibrary::getTet(tetIdx);
//		mTetsIds[tetIdx + 1] = UTetGenFunctionLibrary::getTet(tetIdx + 1);
//		mTetsIds[tetIdx + 2] = UTetGenFunctionLibrary::getTet(tetIdx + 2);
//		mTetsIds[tetIdx + 3] = UTetGenFunctionLibrary::getTet(tetIdx + 3);
//	}
//}

// Called every frame
void AProceduralMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	
	//UFemFunctions::update(mIdAlgoFEM, DeltaTime, &mVerticesBuffer[0]);

	for (int j = 0; j < mVertices.Num(); ++j)
	{
		mVertices[j][0] = mVerticesBuffer[j * 3];
		mVertices[j][1] = mVerticesBuffer[j * 3 + 1];
		mVertices[j][2] = mVerticesBuffer[j * 3 + 2];
	}
	 
	//t = t + GetWorld()->GetDeltaSeconds();
	//	
	//for (int i = 0; i < mVertices.Num(); ++i)
	//{
	//	mVertices[i] = mVertices[i] + FVector(0.2, 0, 0) * t;
	//}
	//
	//UFemFunctions::runFem(mVerticesBuffer, 3 * mVertices.Num(), mTetsBuffer, mTetsIds.Num());
	//
	mMeshComponent->UpdateMeshSection(0, mVertices, mNormals, mUVs, mVertexColors, mTangents);
}
