//==============================================================
// Copyright © 2021 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#pragma once

#ifdef FEMIMP_DLL_EXPORTS
#define FEMIMP_DLL_API __declspec(dllexport)
#else
#define FEMIMP_DLL_API __declspec(dllimport)
#endif

class Fem_Algoritm;

class FEMIMP_DLL_API FEM_Factory {
	
public:
	void create(int ndim, int nodof, int nels);
	void init(float* g_coord, int* g_num, int* in_nf, int in_nn);
	void update();

private: 
	Fem_Algoritm* femAlg;
};

FEMIMP_DLL_API float basicTest(float* verticesBuffer, int bufferSize, int* tetsBuffer, int tetsBufferSize);

FEMIMP_DLL_API void elemStiffnessMatrix(float* verticesBuffer,
	int* tetsBuffer,
	float* loads,
	const int nels,
	const int* loads_nodes_ids,
	const int loadsNodeSize);

FEMIMP_DLL_API void elemStiffnessMatrixReference(float* verticesBuffer, int* tetsBuffer);

FEMIMP_DLL_API void assemblyOfElements(float* verticesBuffer, int nbNodes);

//Funtions to unit test

FEMIMP_DLL_API float invert(float* A, int m);
FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n);