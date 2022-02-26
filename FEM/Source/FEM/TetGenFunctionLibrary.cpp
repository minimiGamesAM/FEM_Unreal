// Fill out your copyright notice in the Description page of Project Settings.


#include "TetGenFunctionLibrary.h"

void* UTetGenFunctionLibrary::DllHandle = nullptr;

bool UTetGenFunctionLibrary::LoadDllHandle()
{
	FString DllFilePath = FPaths::ProjectDir() + "/Plugins/TetraedrosGenerador/TetraGen.dll";
	if (FPaths::FileExists(DllFilePath))
	{
		DllHandle = FPlatformProcess::GetDllHandle(*DllFilePath);
	}

	return DllHandle != nullptr;
}

void UTetGenFunctionLibrary::RunTetGen()
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("runTetGen"));
		if (DllExport)
		{
			typedef void(*runTetGen)();
			runTetGen runTetGenRef = (runTetGen)(DllExport);
			runTetGenRef();
		}
	}
}