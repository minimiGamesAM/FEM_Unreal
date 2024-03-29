//==============================================================
// Copyright � 2021 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#pragma once

#ifdef FEMIMP_DLL_EXPORTS
#define FEMIMP_DLL_API __declspec(dllexport)
#else
#define FEMIMP_DLL_API __declspec(dllimport)
#endif

#include <vector>

template<class T>
class Fem_Algoritm;

template<class T>
class FEM_Factory
{
	
public:
	static int create(int ndim, int nodof, int nels, int nod, int nip, const char* element);
	static void loadedNodes(int id, int* nodes, int loaded_nodes, T* vals);
	static void init(int id, T* g_coord, int* g_num, int* in_nf, int in_nn);
	static void setGravityAcceleration(int id, T gravAcc);
	static void setGravityDirection(int id, T* gravDir);
	static void setDamping(int id, const T fk, const T fm);
	static void setMaterialParams(int id, const T e, const T v, const T gamma);
	static long long update(int id, T dtim, T* verticesBuffer);

private: 
	static std::vector<Fem_Algoritm<T>*> femAlg;
};


template class FEMIMP_DLL_API FEM_Factory<float>;
template class FEMIMP_DLL_API FEM_Factory<double>;

FEMIMP_DLL_API float basicTest(float* verticesBuffer, int bufferSize, int* tetsBuffer, int tetsBufferSize);

FEMIMP_DLL_API void testSparse();

//Funtions to unit test

//FEMIMP_DLL_API float invert(float* A, int m);
//FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n);