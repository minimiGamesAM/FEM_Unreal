// Fill out your copyright notice in the Description page of Project Settings.


#include "ActorTestPlugin.h"

THIRD_PARTY_INCLUDES_START
#include "FemImpLibrary/FemImpLibrary.h"
THIRD_PARTY_INCLUDES_END

// Sets default values
AActorTestPlugin::AActorTestPlugin()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

}

// Called when the game starts or when spawned
void AActorTestPlugin::BeginPlay()
{
	Super::BeginPlay();
	
	float determinante = basicTest();
	GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Determinante %f"), determinante));

}

// Called every frame
void AActorTestPlugin::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

