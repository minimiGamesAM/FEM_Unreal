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

template<class T>
class Fem_Algoritm;

template<class T>
class FEM_Factory
{
	
public:
	void create(int ndim, int nodof, int nels);
	void init(T* g_coord, int* g_num, int* in_nf, int in_nn);
	void update();

private: 
	Fem_Algoritm<T>* femAlg;
};


template class FEMIMP_DLL_API FEM_Factory<float>;
template class FEMIMP_DLL_API FEM_Factory<double>;

FEMIMP_DLL_API float basicTest(float* verticesBuffer, int bufferSize, int* tetsBuffer, int tetsBufferSize);


//Funtions to unit test

//FEMIMP_DLL_API float invert(float* A, int m);
//FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n);