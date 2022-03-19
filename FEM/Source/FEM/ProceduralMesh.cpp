// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralMesh.h"
#include "TetGenFunctionLibrary.h"

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

	UTetGenFunctionLibrary::RunTetGen();
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
		mVertices.Add(p * 100);
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

	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("NB Tets %i "), UTetGenFunctionLibrary::getNumberOfTets()));
	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("NB Facets %i "), UTetGenFunctionLibrary::getNumberOfTrifaces()));
		
	//auto obs = GetDefaultSubobjects<AProceduralMesh>();
	//AProceduralMesh* proMesh = GetComponentByClass<AActor>();
	//
	//if (proMesh)
	//{
	//	proMesh->GetActorGuid();
	//}
}

// Called when the game starts or when spawned
void AProceduralMesh::BeginPlay()
{
	Super::BeginPlay();

	mVertices.Empty();
	mTriangles.Empty();
	mNormals.Empty();
	mUVs.Empty();
	mVertexColors.Empty();
	mTangents.Empty();

	TArray<UActorComponent*> comps;

	this->GetComponents(comps);

	for (int i = 0; i < comps.Num(); ++i) //Because there may be more components
	{
		UProceduralMeshComponent* thisComp = Cast<UProceduralMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisComp)
		{
			mMeshComponent = thisComp;
			auto meshSection = thisComp->GetProcMeshSection(0);

			auto& data = meshSection->ProcVertexBuffer; // .Position;
				
			mVertices.SetNum(data.Num());
			for (int j = 0; j < data.Num(); ++j)
			{
				mVertices[j] = data[j].Position;
			}
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

// Called every frame
void AProceduralMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
	
	t = t + GetWorld()->GetDeltaSeconds();
		
	for (int i = 0; i < mVertices.Num(); ++i)
	{
		mVertices[i] = mVertices[i] + FVector(0.2, 0, 0) * t;
	}

	mMeshComponent->UpdateMeshSection(0, mVertices, mNormals, mUVs, mVertexColors, mTangents);
}
