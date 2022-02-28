#pragma once

#ifdef TETRAGEN_EXPORTS
#define TETRAGEN_API __declspec(dllexport)
#else
#define TETRAGEN_API __declspec(dllimport)
#endif

extern "C" TETRAGEN_API void runTetGen();