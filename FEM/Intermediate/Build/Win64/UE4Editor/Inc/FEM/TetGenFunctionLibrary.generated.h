// Copyright Epic Games, Inc. All Rights Reserved.
/*===========================================================================
	Generated code exported from UnrealHeaderTool.
	DO NOT modify this manually! Edit the corresponding .h files instead!
===========================================================================*/

#include "UObject/ObjectMacros.h"
#include "UObject/ScriptMacros.h"

PRAGMA_DISABLE_DEPRECATION_WARNINGS
#ifdef FEM_TetGenFunctionLibrary_generated_h
#error "TetGenFunctionLibrary.generated.h already included, missing '#pragma once' in TetGenFunctionLibrary.h"
#endif
#define FEM_TetGenFunctionLibrary_generated_h

#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_SPARSE_DATA
#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_RPC_WRAPPERS \
 \
	DECLARE_FUNCTION(execRunTetGen);


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_RPC_WRAPPERS_NO_PURE_DECLS \
 \
	DECLARE_FUNCTION(execRunTetGen);


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_INCLASS_NO_PURE_DECLS \
private: \
	static void StaticRegisterNativesUTetGenFunctionLibrary(); \
	friend struct Z_Construct_UClass_UTetGenFunctionLibrary_Statics; \
public: \
	DECLARE_CLASS(UTetGenFunctionLibrary, UBlueprintFunctionLibrary, COMPILED_IN_FLAGS(0), CASTCLASS_None, TEXT("/Script/FEM"), NO_API) \
	DECLARE_SERIALIZER(UTetGenFunctionLibrary)


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_INCLASS \
private: \
	static void StaticRegisterNativesUTetGenFunctionLibrary(); \
	friend struct Z_Construct_UClass_UTetGenFunctionLibrary_Statics; \
public: \
	DECLARE_CLASS(UTetGenFunctionLibrary, UBlueprintFunctionLibrary, COMPILED_IN_FLAGS(0), CASTCLASS_None, TEXT("/Script/FEM"), NO_API) \
	DECLARE_SERIALIZER(UTetGenFunctionLibrary)


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_STANDARD_CONSTRUCTORS \
	/** Standard constructor, called after all reflected properties have been initialized */ \
	NO_API UTetGenFunctionLibrary(const FObjectInitializer& ObjectInitializer = FObjectInitializer::Get()); \
	DEFINE_DEFAULT_OBJECT_INITIALIZER_CONSTRUCTOR_CALL(UTetGenFunctionLibrary) \
	DECLARE_VTABLE_PTR_HELPER_CTOR(NO_API, UTetGenFunctionLibrary); \
	DEFINE_VTABLE_PTR_HELPER_CTOR_CALLER(UTetGenFunctionLibrary); \
private: \
	/** Private move- and copy-constructors, should never be used */ \
	NO_API UTetGenFunctionLibrary(UTetGenFunctionLibrary&&); \
	NO_API UTetGenFunctionLibrary(const UTetGenFunctionLibrary&); \
public:


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_ENHANCED_CONSTRUCTORS \
	/** Standard constructor, called after all reflected properties have been initialized */ \
	NO_API UTetGenFunctionLibrary(const FObjectInitializer& ObjectInitializer = FObjectInitializer::Get()) : Super(ObjectInitializer) { }; \
private: \
	/** Private move- and copy-constructors, should never be used */ \
	NO_API UTetGenFunctionLibrary(UTetGenFunctionLibrary&&); \
	NO_API UTetGenFunctionLibrary(const UTetGenFunctionLibrary&); \
public: \
	DECLARE_VTABLE_PTR_HELPER_CTOR(NO_API, UTetGenFunctionLibrary); \
	DEFINE_VTABLE_PTR_HELPER_CTOR_CALLER(UTetGenFunctionLibrary); \
	DEFINE_DEFAULT_OBJECT_INITIALIZER_CONSTRUCTOR_CALL(UTetGenFunctionLibrary)


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_PRIVATE_PROPERTY_OFFSET
#define FEM_Source_FEM_TetGenFunctionLibrary_h_12_PROLOG
#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_GENERATED_BODY_LEGACY \
PRAGMA_DISABLE_DEPRECATION_WARNINGS \
public: \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_PRIVATE_PROPERTY_OFFSET \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_SPARSE_DATA \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_RPC_WRAPPERS \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_INCLASS \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_STANDARD_CONSTRUCTORS \
public: \
PRAGMA_ENABLE_DEPRECATION_WARNINGS


#define FEM_Source_FEM_TetGenFunctionLibrary_h_15_GENERATED_BODY \
PRAGMA_DISABLE_DEPRECATION_WARNINGS \
public: \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_PRIVATE_PROPERTY_OFFSET \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_SPARSE_DATA \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_RPC_WRAPPERS_NO_PURE_DECLS \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_INCLASS_NO_PURE_DECLS \
	FEM_Source_FEM_TetGenFunctionLibrary_h_15_ENHANCED_CONSTRUCTORS \
private: \
PRAGMA_ENABLE_DEPRECATION_WARNINGS


template<> FEM_API UClass* StaticClass<class UTetGenFunctionLibrary>();

#undef CURRENT_FILE_ID
#define CURRENT_FILE_ID FEM_Source_FEM_TetGenFunctionLibrary_h


PRAGMA_ENABLE_DEPRECATION_WARNINGS
