// Copyright Epic Games, Inc. All Rights Reserved.
/*===========================================================================
	Generated code exported from UnrealHeaderTool.
	DO NOT modify this manually! Edit the corresponding .h files instead!
===========================================================================*/

#include "UObject/GeneratedCppIncludes.h"
#include "FEM/ProceduralMesh.h"
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable : 4883)
#endif
PRAGMA_DISABLE_DEPRECATION_WARNINGS
void EmptyLinkFunctionForGeneratedCodeProceduralMesh() {}
// Cross Module References
	FEM_API UClass* Z_Construct_UClass_AProceduralMesh_NoRegister();
	FEM_API UClass* Z_Construct_UClass_AProceduralMesh();
	ENGINE_API UClass* Z_Construct_UClass_AActor();
	UPackage* Z_Construct_UPackage__Script_FEM();
// End Cross Module References
	void AProceduralMesh::StaticRegisterNativesAProceduralMesh()
	{
	}
	UClass* Z_Construct_UClass_AProceduralMesh_NoRegister()
	{
		return AProceduralMesh::StaticClass();
	}
	struct Z_Construct_UClass_AProceduralMesh_Statics
	{
		static UObject* (*const DependentSingletons[])();
#if WITH_METADATA
		static const UE4CodeGen_Private::FMetaDataPairParam Class_MetaDataParams[];
#endif
		static const FCppClassTypeInfoStatic StaticCppClassTypeInfo;
		static const UE4CodeGen_Private::FClassParams ClassParams;
	};
	UObject* (*const Z_Construct_UClass_AProceduralMesh_Statics::DependentSingletons[])() = {
		(UObject* (*)())Z_Construct_UClass_AActor,
		(UObject* (*)())Z_Construct_UPackage__Script_FEM,
	};
#if WITH_METADATA
	const UE4CodeGen_Private::FMetaDataPairParam Z_Construct_UClass_AProceduralMesh_Statics::Class_MetaDataParams[] = {
		{ "IncludePath", "ProceduralMesh.h" },
		{ "ModuleRelativePath", "ProceduralMesh.h" },
	};
#endif
	const FCppClassTypeInfoStatic Z_Construct_UClass_AProceduralMesh_Statics::StaticCppClassTypeInfo = {
		TCppClassTypeTraits<AProceduralMesh>::IsAbstract,
	};
	const UE4CodeGen_Private::FClassParams Z_Construct_UClass_AProceduralMesh_Statics::ClassParams = {
		&AProceduralMesh::StaticClass,
		"Engine",
		&StaticCppClassTypeInfo,
		DependentSingletons,
		nullptr,
		nullptr,
		nullptr,
		UE_ARRAY_COUNT(DependentSingletons),
		0,
		0,
		0,
		0x009000A4u,
		METADATA_PARAMS(Z_Construct_UClass_AProceduralMesh_Statics::Class_MetaDataParams, UE_ARRAY_COUNT(Z_Construct_UClass_AProceduralMesh_Statics::Class_MetaDataParams))
	};
	UClass* Z_Construct_UClass_AProceduralMesh()
	{
		static UClass* OuterClass = nullptr;
		if (!OuterClass)
		{
			UE4CodeGen_Private::ConstructUClass(OuterClass, Z_Construct_UClass_AProceduralMesh_Statics::ClassParams);
		}
		return OuterClass;
	}
	IMPLEMENT_CLASS(AProceduralMesh, 2139803760);
	template<> FEM_API UClass* StaticClass<AProceduralMesh>()
	{
		return AProceduralMesh::StaticClass();
	}
	static FCompiledInDefer Z_CompiledInDefer_UClass_AProceduralMesh(Z_Construct_UClass_AProceduralMesh, &AProceduralMesh::StaticClass, TEXT("/Script/FEM"), TEXT("AProceduralMesh"), false, nullptr, nullptr, nullptr);
	DEFINE_VTABLE_PTR_HELPER_CTOR(AProceduralMesh);
PRAGMA_ENABLE_DEPRECATION_WARNINGS
#ifdef _MSC_VER
#pragma warning (pop)
#endif
