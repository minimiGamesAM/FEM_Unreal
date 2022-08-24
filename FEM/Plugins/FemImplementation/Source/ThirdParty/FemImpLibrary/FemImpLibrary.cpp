//==============================================================
// Copyright © 2021 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#include "FemImpLibrary.h"

#include <iostream>

#include "mkl.h"
#include "mkl_pblas.h"
#include <algorithm>

//using namespace sycl;

namespace
{
    int formnf(int* nf, const int nodof, const int nn)
    {
        int m = 0;

        for (int j = 0; j < nn; ++j)
        {
            for (int i = 0; i < nodof; ++i)
            {
                if (nf[j + nn * i] != 0)
                {
                    m = m + 1;
                    nf[j + nn * i] = m;
                }
            }
        }

        return m;
    }
}

float matrixInversion()
{
    lapack_int m = 3;
    float* A = new float[m * m];
    lapack_int* ipiv = new lapack_int[m];

    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;
    A[6] = 7;
    A[7] = 2;
    A[8] = 9;

    LAPACKE_sgetrf(CblasRowMajor, m, m, A, m, ipiv);

    float determinant = 1;

    for (int i = 0; i < m; ++i)
    {
        determinant = determinant * A[i + i * m];
    }

    determinant = -determinant;

    std::cout << "Determinante " << determinant << std::endl;

    LAPACKE_sgetri(CblasRowMajor, m, A, m, ipiv);

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            std::cout << A[j + i * m] << std::endl;
        }
    }
    return determinant;
}

void matrixProductVector()
{
    const MKL_INT   m(3), n(2);
    MKL_INT	lda(n);
    float* a, * x, * y;
    CBLAS_LAYOUT    layout(CblasRowMajor);

    a = new float[m * n]; // (float*)calloc((m) * (n), sizeof(float));
    x = new float[n]; // (float*)calloc(n, sizeof(float));

    y = new float[m];

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            a[j + i * n] = 1.0f;// j + i * n;
        }
    }

    for (int j = 0; j < n; ++j)
    {
        x[j] = 1.0f;
    }

    std::cout << "Matrix antes" << std::endl;

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            std::cout << a[j + i * n] << std::endl;
        }

        y[i] = 0;
    }

    cblas_sgemv(layout, CblasNoTrans, m, n, 1, a, lda, x, 1, 0, y, 1);

    std::cout << "Vector despues" << std::endl;

    for (int j = 0; j < m; ++j)
    {
        std::cout << y[j] << std::endl;
    }

    delete a;
    delete x;
    delete y;

    //PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, a, &lda, "A");
}

float dotProduct(MKL_INT n, float* x, float* y)
{
    return cblas_sdot(n, x, 1, y, 1);
}

float vectorOperation(float* verticesBuffer, int verticesBufferSize)
{
    //MKL_INT  n, incx, incy, i;
    //float* x, * y;
    float    res;
    //MKL_INT  len_x, len_y;
    //
    //n = 1;
    //incx = 1;
    //incy = 1;
    //
    ////len_x = 1 + (n - 1) * abs(incx);
    ////len_y = 1 + (n - 1) * abs(incy);
    //x = new float[n]; // (float*)calloc(len_x, sizeof(float));
    //y = new float[n]; // (float*)calloc(len_y, sizeof(float));
    //if (x == NULL || y == NULL) {
    //    //printf("\n Can't allocate memory for arrays\n");
    //    return;
    //}
    //
    //for (i = 0; i < n; i++) {
    //    x[i] = 2.0f;
    //    y[i] = 5.0f;
    //}

    res = dotProduct(verticesBufferSize, verticesBuffer, verticesBuffer);

    // printf("\n       SDOT = %7.3f", res);
   // std::cout << res << std::endl;
    //free(x);
    //free(y);

    return res;
}

// C = alpha A * B + beta * C
FEMIMP_DLL_API void matmul(float* A, float* B, float* C, int m, int k, int n)
{
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);
}

void matmulTransA(float* A, float* B, float* C, int m, int k, int n)
{
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, B, n, 0, C, n);
}

// mkl_?getrfnp could be an alternative withouh using pivot
FEMIMP_DLL_API float invert(float* A, int m)
{
    lapack_int* ipiv = new lapack_int[m];
        
    auto info = LAPACKE_sgetrf(CblasRowMajor, m, m, A, m, ipiv);
    
    float determinant = 1;
    
    for (int i = 0; i < m; ++i)
    {
        determinant = determinant * A[i + i * m];
        determinant = ipiv[i] - 1 != i ? -determinant : determinant;
    }
            
    LAPACKE_sgetri(CblasRowMajor, m, A, m, ipiv);
    
    return determinant;
}

void beemat(float* bee, int nbRowsBee, int nbColumsBee, float* deriv, int nbColumsDerivs)
{
    int k, l, n, ih(nbRowsBee), nod(nbColumsDerivs);
    float x, y, z;

    //ih = UBOUND(bee, 1) 
    //nod = UBOUND(deriv, 2)

    for (int m = 1; m < nod + 1; ++m)
    {
        n = 3 * m;
        k = n - 1;
        l = k - 1;

        x = deriv[m - 1 + 0 * nbColumsDerivs]; // (1, m);
        y = deriv[m - 1 + 1 * nbColumsDerivs]; //(2, m);
        z = deriv[m - 1 + 2 * nbColumsDerivs]; //(3, m);

        bee[l - 1 + 0 * nbColumsBee] = x; //(1, l) = x
        bee[k - 1 + 3 * nbColumsBee] = x; //(4, k) = x
        bee[n - 1 + 5 * nbColumsBee] = x; //(6, n) = x
        bee[k - 1 + 1 * nbColumsBee] = y; //(2, k) = y
        bee[l - 1 + 3 * nbColumsBee] = y; //(4, l) = y
        bee[n - 1 + 4 * nbColumsBee] = y; //(5, n) = y
        bee[n - 1 + 2 * nbColumsBee] = z; //(3, n) = z
        bee[k - 1 + 4 * nbColumsBee] = z; //(5, k) = z
        bee[l - 1 + 5 * nbColumsBee] = z; //(6, l) = z
    }
}

void deemat(float* dee, int nbRowsDee, float e, float v)
{
    float ih = nbRowsDee;
    float one(1.0f), two(2.0f), pt5 = 0.5f;
    float v1 = one - v;
    float c = e / ((one + v) * (one - two * v));
    
    float v2 = v / (one - v);
    float vv = (one - two * v) / (one - v) * pt5;

    float factor = e / (two * (one + v) * vv);

    for (int i = 1; i <= 3; ++i)
    {
        dee[i - 1 + (i - 1) * nbRowsDee] = one * factor;
    }

    for (int i = 4; i <= 6; ++i)
    {
        dee[i - 1 + (i - 1) * nbRowsDee] = vv * factor;
    }

    v2 *= factor;
    dee[1 + 0 * nbRowsDee] = v2; //dee[1, 2] = v2  
    dee[0 + 1 * nbRowsDee] = v2; //dee[2, 1] = v2
    dee[2 + 0 * nbRowsDee] = v2; //dee[1, 3] = v2
    dee[0 + 2 * nbRowsDee] = v2; //dee[3, 1] = v2
    dee[2 + 1 * nbRowsDee] = v2; //dee[2, 3] = v2
    dee[1 + 2 * nbRowsDee] = v2; //dee[3, 2] = v2
}                                         

FEMIMP_DLL_API void elemStiffnessMatrix(float* g_coord, int* g_num)
{
    //nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
    const int nodof = 3;
    //nn = total number of nodes in the problem
    const int nn = 8;

    //nf = nodal freedom array(nodof rows and nn colums)
    int nf[nodof * nn] = {  0, 1, 0, 1, 0, 0, 1, 1,
                            0, 0, 0, 0, 1, 1, 1, 1,
                            1, 1, 0, 0, 1, 0, 0, 1 };
      
    //neq = number of degree of freedom in the mesh
    int neq = formnf(nf, nodof, nn);

    for (int i = 0; i < nodof * nn; ++i)
    {
    	std::cout << "nf " << nf[i] << std::endl;
    }

    std::cout << neq << std::endl;
    //nod = number of node per element
    const int nod = 4;
    //ndof = number of degree of freedom per element
    const int ndof = nod * nodof;


}

FEMIMP_DLL_API void elemStiffnessMatrixReference(float* verticesBuffer, int* tetsBuffer)
{
    const int nbPoints = 4;
    const int dim = 3;

    const float e = 10000.0f;
    const float v = 0.3f;
    
    float* coord = new float[nbPoints * dim];

    for (int tetIdCoord = 0; tetIdCoord < nbPoints; ++tetIdCoord)
    {       
        int pointId = tetsBuffer[tetIdCoord];

        for (int j = 0; j < dim; ++j)
        {
            coord[j + dim * tetIdCoord] = verticesBuffer[dim * pointId + j];
        }
    }

    // derivate of the shape function
    float der[dim * nbPoints] = {
        1.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, -1.0f,
        0.0f, 0.0f, 1.0f, -1.0f
    };

    // calculate jac
    float jac[dim * dim] = {};
    matmul(der, coord, jac, dim, nbPoints, dim);
    float det = invert(jac, dim);

    // calculate the derivative in x, y, z
    float deriv[dim * nbPoints] = {};
       
    matmul(jac, der, deriv, dim, dim, nbPoints);

    const int beeBufferSize = 6 * (nbPoints * 3);
    const int deeBufferSize = 6 * 6;
    //float* bee = new float[beeBufferSize];// = {};
    
    float bee[beeBufferSize] = {};
    float* beeItt = bee;
    std::fill(beeItt, beeItt + beeBufferSize, 0.0f);

    float dee[deeBufferSize] = {};
    float* deeItt = dee;
    std::fill(deeItt, deeItt + deeBufferSize, 0.0f);

    beemat(bee, 6, (nbPoints * 3), deriv, nbPoints);
    deemat(dee, 6, e, v);

    float temporal[beeBufferSize] = {};
    float* tempItt = temporal;
    std::fill(tempItt, tempItt + beeBufferSize, 0.0f);

    const int btdbBufferSize = (nbPoints * 3) * (nbPoints * 3);
    float btdb[btdbBufferSize] = {};
    float* btdbItt = btdb;
    std::fill(btdbItt, btdbItt + btdbBufferSize, 0.0f);

    matmulTransA(bee, dee, temporal, nbPoints * 3, 6, 6);
    matmul(temporal, bee, btdb, nbPoints * 3, 6, nbPoints * 3);
    
    for (int i = 0; i < btdbBufferSize; ++i)
    {
        btdb[i] *= det / 6;
    }

    for (int i = 0; i < btdbBufferSize; ++i)
    {
        std::cout << btdb[i] << std::endl;
    }
}

FEMIMP_DLL_API void assemblyOfElements(float* verticesBuffer, const int nbNodes)
{
   
}


FEMIMP_DLL_API float basicTest(float* verticesBuffer, int verticesBufferSize, int* tetsBuffer, int tetsBufferSize)
{
    float determinat = matrixInversion();
    matrixProductVector();
    float dotResult = vectorOperation(verticesBuffer, verticesBufferSize);

    std::cout << std::endl;
    std::cout << "GIRANDO TEST BASICO" << std::endl;
    //MessageBox(NULL, TEXT("BASIC DLL CARGADO INTEL."), TEXT("Third Party Plugin"), MB_OK);
    return dotResult;
}