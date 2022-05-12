// Copyright Epic Games, Inc. All Rights Reserved.

#include "FemImplementation.h"
#include "Core.h"
#include "Interfaces/IPluginManager.h"
#include "FemImpLibrary/FemImpLibrary.h"

#define LOCTEXT_NAMESPACE "FFemImplementationModule"

void FFemImplementationModule::StartupModule()
{
	// This code will execute after your module is loaded into memory; the exact timing is specified in the .uplugin file per-module
	
	FString BaseDir = IPluginManager::Get().FindPlugin("FemImplementation")->GetBaseDir();
	FString LibraryPath;


	LibraryPath = FPaths::Combine(*BaseDir, TEXT("Binaries/ThirdParty/FemImpLibrary/Win64/mkl_avx2.2.dll"));
	ExampleLibraryHandle = !LibraryPath.IsEmpty() ? FPlatformProcess::GetDllHandle(*LibraryPath) : nullptr;

	LibraryPath = FPaths::Combine(*BaseDir, TEXT("Binaries/ThirdParty/FemImpLibrary/Win64/mkl_core.2.dll"));
	ExampleLibraryHandle = !LibraryPath.IsEmpty() ? FPlatformProcess::GetDllHandle(*LibraryPath) : nullptr;
	LibraryPath = FPaths::Combine(*BaseDir, TEXT("Binaries/ThirdParty/FemImpLibrary/Win64/mkl_sequential.2.dll"));
	ExampleLibraryHandle = !LibraryPath.IsEmpty() ? FPlatformProcess::GetDllHandle(*LibraryPath) : nullptr;
	
	LibraryPath = FPaths::Combine(*BaseDir, TEXT("Binaries/ThirdParty/FemImpLibrary/Win64/FemImpLibrary.dll"));

	ExampleLibraryHandle = !LibraryPath.IsEmpty() ? FPlatformProcess::GetDllHandle(*LibraryPath) : nullptr;

	if (ExampleLibraryHandle)
	{

		float determinante = basicTest();
		auto message = FString("ThirdPartyLibraryError %f", determinante);// +FString(determinante);
		//FMessageDialog::Open(EAppMsgType::Ok, FText(message));
		//GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Yellow, FString::Printf(TEXT("Determinante %f"), determinante));
	}
}

void FFemImplementationModule::ShutdownModule()
{
	// This function may be called during shutdown to clean up your module.  For modules that support dynamic reloading,
	// we call this function before unloading the module.
}

#undef LOCTEXT_NAMESPACE
	
IMPLEMENT_MODULE(FFemImplementationModule, FemImplementation)