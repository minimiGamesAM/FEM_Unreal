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

FEMIMP_DLL_API float basicTest(float* verticesBuffer, int bufferSize, int* tetsBuffer, int tetsBufferSize);

FEMIMP_DLL_API void elemStiffnessMatrix(float* verticesBuffer, int* tetsBuffer);

FEMIMP_DLL_API void elemStiffnessMatrixReference(float* verticesBuffer, int* tetsBuffer);

FEMIMP_DLL_API void assemblyOfElements(float* verticesBuffer, int nbNodes);

//Funtions to unit test

FEMIMP_DLL_API float invert(float* A, int m);
FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n);