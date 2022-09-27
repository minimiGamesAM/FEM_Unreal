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

#include <vector>

template<class T>
class Fem_Algoritm;

template<class T>
class FEM_Factory
{
	
public:
	static int create(int ndim, int nodof, int nels);
	static void init(int id, T* g_coord, int* g_num, int* in_nf, int in_nn);
	static void update(int id, T dtim, T* verticesBuffer);

private: 
	static std::vector<Fem_Algoritm<T>*> femAlg;
};


template class FEMIMP_DLL_API FEM_Factory<float>;
template class FEMIMP_DLL_API FEM_Factory<double>;

FEMIMP_DLL_API float basicTest(float* verticesBuffer, int bufferSize, int* tetsBuffer, int tetsBufferSize);


//Funtions to unit test

//FEMIMP_DLL_API float invert(float* A, int m);
//FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n);