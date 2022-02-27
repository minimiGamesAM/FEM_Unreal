// Copyright Epic Games, Inc. All Rights Reserved.
/*===========================================================================
	Generated code exported from UnrealHeaderTool.
	DO NOT modify this manually! Edit the corresponding .h files instead!
===========================================================================*/

#include "UObject/GeneratedCppIncludes.h"
#include "FEM/TetGenFunctionLibrary.h"
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable : 4883)
#endif
PRAGMA_DISABLE_DEPRECATION_WARNINGS
void EmptyLinkFunctionForGeneratedCodeTetGenFunctionLibrary() {}
// Cross Module References
	FEM_API UClass* Z_Construct_UClass_UTetGenFunctionLibrary_NoRegister();
	FEM_API UClass* Z_Construct_UClass_UTetGenFunctionLibrary();
	ENGINE_API UClass* Z_Construct_UClass_UBlueprintFunctionLibrary();
	UPackage* Z_Construct_UPackage__Script_FEM();
// End Cross Module References
	DEFINE_FUNCTION(UTetGenFunctionLibrary::execRunTetGen)
	{
		P_FINISH;
		P_NATIVE_BEGIN;
		UTetGenFunctionLibrary::RunTetGen();
		P_NATIVE_END;
	}
	void UTetGenFunctionLibrary::StaticRegisterNativesUTetGenFunctionLibrary()
	{
		UClass* Class = UTetGenFunctionLibrary::StaticClass();
		static const FNameNativePtrPair Funcs[] = {
			{ "RunTetGen", &UTetGenFunctionLibrary::execRunTetGen },
		};
		FNativeFunctionRegistrar::RegisterFunctions(Class, Funcs, UE_ARRAY_COUNT(Funcs));
	}
	struct Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics
	{
#if WITH_METADATA
		static const UE4CodeGen_Private::FMetaDataPairParam Function_MetaDataParams[];
#endif
		static const UE4CodeGen_Private::FFunctionParams FuncParams;
	};
#if WITH_METADATA
	const UE4CodeGen_Private::FMetaDataPairParam Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics::Function_MetaDataParams[] = {
		{ "ModuleRelativePath", "TetGenFunctionLibrary.h" },
	};
#endif
	const UE4CodeGen_Private::FFunctionParams Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics::FuncParams = { (UObject*(*)())Z_Construct_UClass_UTetGenFunctionLibrary, nullptr, "RunTetGen", nullptr, nullptr, 0, nullptr, 0, RF_Public|RF_Transient|RF_MarkAsNative, (EFunctionFlags)0x04022401, 0, 0, METADATA_PARAMS(Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics::Function_MetaDataParams, UE_ARRAY_COUNT(Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics::Function_MetaDataParams)) };
	UFunction* Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen()
	{
		static UFunction* ReturnFunction = nullptr;
		if (!ReturnFunction)
		{
			UE4CodeGen_Private::ConstructUFunction(ReturnFunction, Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen_Statics::FuncParams);
		}
		return ReturnFunction;
	}
	UClass* Z_Construct_UClass_UTetGenFunctionLibrary_NoRegister()
	{
		return UTetGenFunctionLibrary::StaticClass();
	}
	struct Z_Construct_UClass_UTetGenFunctionLibrary_Statics
	{
		static UObject* (*const DependentSingletons[])();
		static const FClassFunctionLinkInfo FuncInfo[];
#if WITH_METADATA
		static const UE4CodeGen_Private::FMetaDataPairParam Class_MetaDataParams[];
#endif
		static const FCppClassTypeInfoStatic StaticCppClassTypeInfo;
		static const UE4CodeGen_Private::FClassParams ClassParams;
	};
	UObject* (*const Z_Construct_UClass_UTetGenFunctionLibrary_Statics::DependentSingletons[])() = {
		(UObject* (*)())Z_Construct_UClass_UBlueprintFunctionLibrary,
		(UObject* (*)())Z_Construct_UPackage__Script_FEM,
	};
	const FClassFunctionLinkInfo Z_Construct_UClass_UTetGenFunctionLibrary_Statics::FuncInfo[] = {
		{ &Z_Construct_UFunction_UTetGenFunctionLibrary_RunTetGen, "RunTetGen" }, // 1495578352
	};
#if WITH_METADATA
	const UE4CodeGen_Private::FMetaDataPairParam Z_Construct_UClass_UTetGenFunctionLibrary_Statics::Class_MetaDataParams[] = {
		{ "Comment", "/**\n * \n */" },
		{ "IncludePath", "TetGenFunctionLibrary.h" },
		{ "ModuleRelativePath", "TetGenFunctionLibrary.h" },
	};
#endif
	const FCppClassTypeInfoStatic Z_Construct_UClass_UTetGenFunctionLibrary_Statics::StaticCppClassTypeInfo = {
		TCppClassTypeTraits<UTetGenFunctionLibrary>::IsAbstract,
	};
	const UE4CodeGen_Private::FClassParams Z_Construct_UClass_UTetGenFunctionLibrary_Statics::ClassParams = {
		&UTetGenFunctionLibrary::StaticClass,
		nullptr,
		&StaticCppClassTypeInfo,
		DependentSingletons,
		FuncInfo,
		nullptr,
		nullptr,
		UE_ARRAY_COUNT(DependentSingletons),
		UE_ARRAY_COUNT(FuncInfo),
		0,
		0,
		0x001000A0u,
		METADATA_PARAMS(Z_Construct_UClass_UTetGenFunctionLibrary_Statics::Class_MetaDataParams, UE_ARRAY_COUNT(Z_Construct_UClass_UTetGenFunctionLibrary_Statics::Class_MetaDataParams))
	};
	UClass* Z_Construct_UClass_UTetGenFunctionLibrary()
	{
		static UClass* OuterClass = nullptr;
		if (!OuterClass)
		{
			UE4CodeGen_Private::ConstructUClass(OuterClass, Z_Construct_UClass_UTetGenFunctionLibrary_Statics::ClassParams);
		}
		return OuterClass;
	}
	IMPLEMENT_CLASS(UTetGenFunctionLibrary, 1308876088);
	template<> FEM_API UClass* StaticClass<UTetGenFunctionLibrary>()
	{
		return UTetGenFunctionLibrary::StaticClass();
	}
	static FCompiledInDefer Z_CompiledInDefer_UClass_UTetGenFunctionLibrary(Z_Construct_UClass_UTetGenFunctionLibrary, &UTetGenFunctionLibrary::StaticClass, TEXT("/Script/FEM"), TEXT("UTetGenFunctionLibrary"), false, nullptr, nullptr, nullptr);
	DEFINE_VTABLE_PTR_HELPER_CTOR(UTetGenFunctionLibrary);
PRAGMA_ENABLE_DEPRECATION_WARNINGS
#ifdef _MSC_VER
#pragma warning (pop)
#endif
