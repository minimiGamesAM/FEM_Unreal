//==============================================================
// Copyright © 2021 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#pragma once

#ifdef SENCILLODLL_EXPORTS
#define DPCPP_DLL_API __declspec(dllexport)
#else
#define DPCPP_DLL_API __declspec(dllimport)
#endif

#include <array>
#include <iostream>

// ARRAY type & data size for use in this example
//constexpr size_t array_size = 10000;
//using IntArray=std::array<int, array_size>;

DPCPP_DLL_API void matrixInversion();
DPCPP_DLL_API void matrixProductVector();
DPCPP_DLL_API void vectorOperation();
