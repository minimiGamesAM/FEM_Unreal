#pragma once

#ifdef TETRAGEN_EXPORTS
#define TETRAGEN_API __declspec(dllexport)
#else
#define TETRAGEN_API __declspec(dllimport)
#endif

extern "C" TETRAGEN_API void runTetGen(char* file, char* s);

extern "C" TETRAGEN_API void runTetGen2(float* points, int* faces, int facesSize, int pointsSizes, char* switches);

extern "C" TETRAGEN_API int getNumberOfPoints();

extern "C" TETRAGEN_API double getPoint(int idx);

extern "C" TETRAGEN_API int getNumberOfTrifaces();

extern "C" TETRAGEN_API int getTrifacet(int idx);

extern "C" TETRAGEN_API int getNumberOfTets();

extern "C" TETRAGEN_API int getTet(int idx);

extern "C" TETRAGEN_API int getTet2facelist(int tetIdx, int idx);