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

//https://www.youtube.com/watch?v=dKlMEmVgbvg

void UTetGenFunctionLibrary::RunTetGen()
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("runTetGen"));
		if (DllExport)
		{
			typedef void(*runTetGen)(char* file_poly, char* switches);
			runTetGen runTetGenRef = (runTetGen)(DllExport);
			
			FString fileModel = FPaths::ProjectDir() + "/Content/Models/verificacion";
			
			char* result = TCHAR_TO_ANSI(*fileModel);
						
			runTetGenRef(result, "pqz-f-nn");// 
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

float UTetGenFunctionLibrary::getPoint(int idx)
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getPoint"));
		if (DllExport)
		{
			typedef double(*getPoint)(int idx);
			getPoint getPointRef = (getPoint)(DllExport);

			return float(getPointRef(idx));
		}
	}

	return 0.0f;
}

int UTetGenFunctionLibrary::getNumberOfTrifaces()
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getNumberOfTrifaces"));
		if (DllExport)
		{
			typedef int(*getNumberOfTrifaces)();
			getNumberOfTrifaces getNumberOfTrifacesRef = (getNumberOfTrifaces)(DllExport);

			return getNumberOfTrifacesRef();
		}
	}

	return 0;
}

int UTetGenFunctionLibrary::getTrifacet(int idx)
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getTrifacet"));
		if (DllExport)
		{
			typedef int(*getTrifacet)(int idx);
			getTrifacet getTrifacetRef = (getTrifacet)(DllExport);

			return getTrifacetRef(idx);
		}
	}

	return 0;
}

int UTetGenFunctionLibrary::getNumberOfTets()
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getNumberOfTets"));
		if (DllExport)
		{
			typedef int(*getNumberOfTets)();
			getNumberOfTets getNumberOfTetsRef = (getNumberOfTets)(DllExport);

			return getNumberOfTetsRef();
		}
	}

	return 0;
}

int UTetGenFunctionLibrary::getTet(int idx)
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getTet"));
		if (DllExport)
		{
			typedef int(*getTet)(int idx);
			getTet getTetRef = (getTet)(DllExport);

			return getTetRef(idx);
		}
	}

	return 0;
}

int UTetGenFunctionLibrary::getTet2facelist(int tetIdx, int idx)
{
	if (DllHandle || LoadDllHandle()) //We have a valid dll handle
	{
		void* DllExport = FPlatformProcess::GetDllExport(DllHandle, *FString("getTet2facelist"));
		if (DllExport)
		{
			typedef int(*getTet2facelist)(int tetIdx, int idx);
			getTet2facelist getTet2facelistRef = (getTet2facelist)(DllExport);

			return getTet2facelistRef(tetIdx, idx);
		}
	}

	return 0;
}