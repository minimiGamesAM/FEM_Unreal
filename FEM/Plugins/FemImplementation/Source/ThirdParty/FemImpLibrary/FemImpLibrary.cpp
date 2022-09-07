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
#include <functional>
#include <iterator>

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

    void num_to_g(const int* num, const int* nf, int* g, const int nod, const int nodof, const int nn)
    {
        for (int j = 0; j < nod; ++j)
        {
            for (int k = 0; k < nodof; ++k)
            {
                g[k + j * nodof] = nf[num[j] + k * nn];
            }
        }
    }

    //This subroutine computes the skyline profile.
    void fkdiag(int* kdiag, const int* g, const int ndof)
    {
        int idof = ndof;
        for (int i = 0; i < idof; ++i)
        {
            int iwp1 = 1;

            if (g[i] != 0)
            {
                for (int j = 0; j < idof; ++j)
                {
                    if (g[j] != 0)
                    {
                        int im = g[i] - g[j] + 1;
                        if (im > iwp1) iwp1 = im;
                    }
                }

                int k = g[i] - 1;
                if (iwp1 > kdiag[k]) kdiag[k] = iwp1;
            }
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
        
        for (int i = 1; i <= 3; ++i)
        {
            dee[i - 1 + (i - 1) * nbRowsDee] = one;
        }

        for (int i = 4; i <= 6; ++i)
        {
            dee[i - 1 + (i - 1) * nbRowsDee] = vv;
        }

        //v2 *= factor;
        dee[1 + 0 * nbRowsDee] = v2; //dee[1, 2] = v2  
        dee[0 + 1 * nbRowsDee] = v2; //dee[2, 1] = v2
        dee[2 + 0 * nbRowsDee] = v2; //dee[1, 3] = v2
        dee[0 + 2 * nbRowsDee] = v2; //dee[3, 1] = v2
        dee[2 + 1 * nbRowsDee] = v2; //dee[2, 3] = v2
        dee[1 + 2 * nbRowsDee] = v2; //dee[3, 2] = v2

        for (int i = 0; i < nbRowsDee* nbRowsDee; ++i)
        {
            dee[i] *= e / (two * (one + v) * vv);
        }
    }

    void fsparv(float* kv, const float* km, int* g, const int* kdiag, const int idof)
    {
        for (int i = 0; i < idof; ++i)
        {
            int k = g[i];

            if(k != 0)
            {
                for(int j = 0; j < idof; ++j)
                {
                    if(g[j] != 0)
                    {
                        int iw = k - g[j];
                        if (iw >= 0)
                        {
                            int ival = kdiag[k - 1] - iw - 1;
                            kv[ival] = kv[ival] + km[i + j * idof];
                        }
                    }
                }
            }

        }
    }

    void sparin(float* kv, int* kdiag, int neq)
    {
        int n = neq;
        kv[0] = std::sqrt(kv[0]);

        for (int i = 2; i <= n; ++i)
        {
            float x = 0.0f;
            int ki = kdiag[i - 1] - i;
            int l = kdiag[i - 1 - 1] - ki + 1;

            for (int j = l; j <= i; ++j)
            {
                x = kv[ki + j - 1];
                int kj = kdiag[j - 1] - j;

                if (j != 1)
                {
                    int ll = kdiag[j - 1 - 1] - kj + 1;
                    ll = std::max(l, ll);

                    if (ll != j)
                    {
                        int m = j - 1;

                        for (int k = ll; k <= m; ++k)
                        {
                            x = x - kv[ki + k - 1] * kv[kj + k - 1];
                        }
                    }
                }

                kv[ki + j - 1] = x / kv[kj + j - 1];
            }

            kv[ki + i - 1] = std::sqrt(x);
        }
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

FEMIMP_DLL_API void elemStiffnessMatrix(float* g_coord, int* g_num, float* loads, const int nels)
{
    //nodof = number of freedoms per node (x, y, z, q1, q2, q3 etc)
    const int nodof = 3;
    //nn = total number of nodes in the problem
    const int nn = 8;

    //nip = number of intregation points per element
    const int nip = 1;

    //ndim = number of dimensions
    const int ndim = 3;

    //nprops = number of material properties
    const int nprops = 3;

    //np_types = number of diffent property types
    const int np_types = 1;

    //nf = nodal freedom array(nodof rows and nn colums)
    int nf[nodof * nn] = { 0, 1, 0, 1, 0, 1, 0, 1,
                            0, 0, 0, 0, 1, 1, 1, 1,
                            1, 1, 0, 0, 1, 1, 0, 0 };

    //neq = number of degree of freedom in the mesh
    int neq = formnf(nf, nodof, nn);

    //nod = number of node per element
    const int nod = 4;
    //ndof = number of degree of freedom per element
    const int ndof = nod * nodof;

    //prop = material property(e, v, gamma)
    float* prop = new float[nprops * np_types];
    
    prop[0] = 100.0f;
    prop[1] = 0.3f;
    prop[2] = 0.0f;

    int* g_g = new int[ndof * nels];

    // skyline profile
    int* kdiag = new int[neq];
    std::fill(kdiag, kdiag + neq, 0);

    //---------------------- loop the elements to find global arrays sizes---- -

    for (int i = 0; i < nels; ++i)
    {
        int g[ndof] = {};
        int* num = &g_num[4 * i];

        num_to_g(num, nf, g, nod, nodof, nn);

        for (int j = 0; j < ndof; ++j)
        {
            g_g[i + j * nels] = g[j];
        }

        fkdiag(kdiag, g, ndof);
    }

    for (int i = 1; i < neq; ++i)
    {
        kdiag[i] = kdiag[i] + kdiag[i - 1];
    }

    //kv = global stiffness matrix
    float* kv = new float[kdiag[neq - 1]];
    std::fill(kv, kv + kdiag[neq - 1], 0.0f);

    //gravlo = global gravity loading vector 
    float* gravlo = new float[neq];
    std::fill(gravlo, gravlo + neq, 0.0f);

    //----------------------- element stiffness integration and assembly--------


    //call sample, but for tet there is just one point
    float points[] = { 0.25f, 0.25f, 0.25f };
    float weights = 1.0f / 6.0f;

    //nst = number of stress / strain terms
    const int nst = 6;
    const float e = 100.0f;
    const float v = 0.3f;

    int* etype = new int[nels];
    std::fill(etype, etype + nels, 1);

    //km = element stiffness matrix
    float km[ndof * ndof] = {};

    //eld = element displacement vector
    float eld[ndof] = {};

    float fun[4] = { 0.25f, 0.25f, 0.25f, 0.25f };
        
    float der[ndim * nod] = {
            1.0f, 0.0f, 0.0f, -1.0f,
            0.0f, 1.0f, 0.0f, -1.0f,
            0.0f, 0.0f, 1.0f, -1.0f
    };

    float jac[ndim * ndim] = {};
    float deriv[ndim * nod] = {};
    float bee[nst * ndof] = {};

    for (int i = 0; i < nels; ++i)
    {
        float dee[nst * nst] = {};
        std::fill(std::begin(dee), std::end(dee), 0.0f);
        deemat(dee, nst, e, v);

        int* num = &g_num[4 * i];

        float coord[nod * ndim] = {};

        for (int j = 0; j < nod; ++j)
        {
            for (int k = 0; k < ndim; ++k)
            {
                coord[k + j * ndim] = g_coord[num[j] + k * nn];
            }
        }

        int g[ndof] = {};
        for (int j = 0; j < ndof; ++j)
        {
            g[j] = g_g[i + j * nels];
        }

        std::fill(std::begin(km), std::end(km), 0.0f);
        std::fill(std::begin(eld), std::end(eld), 0.0f);

        ////// for each point of integration, we have just one for 3d tets
        // calculate jac
        matmul(der, coord, jac, ndim, nod, ndim);
        float det = invert(jac, ndim);

        // calculate the derivative in x, y, z
        matmul(jac, der, deriv, ndim, ndim, nod);

        std::fill(std::begin(bee), std::end(bee), 0.0f);
        beemat(bee, nst, ndof, deriv, nod);

        float temporal[nst * ndof] = {};
        std::fill(std::begin(temporal), std::end(temporal), 0.0f);

        matmulTransA(bee, dee, temporal, ndof, nst, nst);
        matmul(temporal, bee, km, ndof, nst, ndof);

        std::for_each(std::begin(km), std::end(km), [&](float& v) { v *= det * weights; });

        for (int j = nodof; j <= ndof; j += nodof)
        {
            eld[j - 1] = fun[int(j / nodof) - 1] * det * weights;
        }

        ////// end for each point

        fsparv(kv, km, g, kdiag, ndof);

        for (int j = 0; j < nels; ++j)
        {
            gravlo[j] += -eld[j] * prop[etype[i] - 1 + np_types * 2];
        }
    }

    // TOBE implemented
    //IF(fixed_freedoms /= 0)THEN
    //    ALLOCATE(node(fixed_freedoms), sense(fixed_freedoms), &
    //        value(fixed_freedoms), no(fixed_freedoms))
    //    READ(10, *)(node(i), sense(i), value(i), i = 1, fixed_freedoms)
    //    DO  i = 1, fixed_freedoms
    //    no(i) = nf(sense(i), node(i))
    //    END DO
    //    kv(kdiag(no)) = kv(kdiag(no)) + penalty
    //    loads(no) = kv(kdiag(no)) * value
    //END IF

    sparin(kv, kdiag, neq);

    //std::for_each(kv, kv + neq, [](float v) {
    //
    //    std::cout << "kv sparin " << v << std::endl;
    //    });

    //for (int j = 0; j < 69; ++j)
    //{
    //    //nprops * np_types
    //    std::cout << "kv ooo " << kv[j] << std::endl;
    //}


    //loads = global load (displacement) vector
    //float* loads = new float[neq];
    //std::fill(loads, loads, 0.0f);


        //for (int j = 0; j < nels; ++j)
        //{
        //    //nprops * np_types
        //    std::cout << "gravlo " << gravlo[j] << std::endl;
        //}

        //std::for_each(kv, kv + kdiag[neq - 1], [](float v) {
        //
        //    std::cout << "kv " << v << std::endl;
        //    });
        //
        //std::cout << "****" << std::endl;
    

    //for (int i = 0; i < neq; ++i)
    //{
    //    std::cout << "kdiag " << i << " " << kdiag[i] << std::endl;
    //}
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