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
			typedef void(*runTetGen)(char* file_poly, char* switches);
			runTetGen runTetGenRef = (runTetGen)(DllExport);
			
			FString fileModel = FPaths::ProjectDir() + "/Content/Models/triangulatedModel";
			
			char* result = TCHAR_TO_ANSI(*fileModel);
						
			runTetGenRef(result, "pq1.414a0.1");
		}
	}
}

int UTetGenFunctionLibrary::getNumberOfPoints()
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getNumberOfPoints"));
		if (DllExport)
		{
			typedef int(*getNumberOfPoints)();
			getNumberOfPoints getNumberOfPointsRef = (getNumberOfPoints)(DllExport);
			
			return getNumberOfPointsRef();
		}
	}

	return 0;
}