#pragma once

#ifdef TETRAGEN_EXPORTS
#define TETRAGEN_API __declspec(dllexport)
#else
#define TETRAGEN_API __declspec(dllimport)
#endif

extern "C" TETRAGEN_API void runTetGen(const float* points,
										const int* faces,
										const int facesSize,
										const int pointsSizes,
										int&     numberOfPoints,
										double*& pointlist,
										int&     numberoftrifaces,
										int*&    trifacelist,
										int&     numberoftetrahedra,
										int*&    tetrahedronlist,
										int*&    tet2facelist,
										char*   switches);
