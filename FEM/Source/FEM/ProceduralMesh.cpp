// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralMesh.h"

//https://nerivec.github.io/old-ue4-wiki/pages/procedural-mesh-component-in-cgetting-started.html
// without the "Pluging part"
// 1) Add in MyProject.Build.cs "ProceduralMeshComponent" at the end of PublicDependencyModuleNames.AddRange
// 2) Add in AdditionalDependencies in .uproject file "ProceduralMeshComponent"
// 3) To fix errors with Visual Studio IntelliSense you need to right-click MyProject.uproject and re-generate Visual Studio project files. In Visual Studio 2017, open "Solution Explorer" and open the "Game" folder, right-click on the first line, which should be the root of your solution, select: "Rescan Solution".
#include "ProceduralMeshComponent.h"

// Sets default values
AProceduralMesh::AProceduralMesh()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;
	
}

// Called when the game starts or when spawned
void AProceduralMesh::BeginPlay()
{
	Super::BeginPlay();

	TArray<UActorComponent*> comps;

	this->GetComponents(comps);

	for (int i = 0; i < comps.Num(); ++i) //Because there may be more components
	{
		UProceduralMeshComponent* thisComp = Cast<UProceduralMeshComponent>(comps[i]); //try to cast to static mesh component
		if (thisComp)
		{
			thisComp->ClearAllMeshSections();
			//This is the static mesh component
			//auto c = thisComp->GetActorGuid();
		
		}
	}

	//auto obs = GetDefaultSubobjects<AProceduralMesh>();
	//AProceduralMesh* proMesh = GetComponentByClass<AActor>();
	//
	//if (proMesh)
	//{
	//	proMesh->GetActorGuid();
	//}

}

// Called every frame
void AProceduralMesh::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

