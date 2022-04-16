//==============================================================
// Copyright © Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
//#include <CL/sycl.hpp>
#include <vector>
#include <iostream>


#include "mkl.h"
#include "mkl_pblas.h"

//using namespace sycl;

void matrixInversion()
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
}

void matrixProductVector()
{
    const MKL_INT   m(3), n(2);
    MKL_INT	lda(n);
    float* a, *x, *y;
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

int vectorOperation()
{
    MKL_INT  n, incx, incy, i;
    float* x, * y;
    float    res;
    MKL_INT  len_x, len_y;

    n = 1;
    incx = 1;
    incy = 1;

    len_x = 1 + (n - 1) * abs(incx);
    len_y = 1 + (n - 1) * abs(incy);
    x = (float*)calloc(len_x, sizeof(float));
    y = (float*)calloc(len_y, sizeof(float));
    if (x == NULL || y == NULL) {
        printf("\n Can't allocate memory for arrays\n");
        return 1;
    }

    for (i = 0; i < n; i++) {
        x[i] = 2.0f;
        y[i] = 1.0f;
    }

    res = dotProduct(n, x, y);

    printf("\n       SDOT = %7.3f", res);

    free(x);
    free(y);

    return 0;
}


//************************************
// Demonstrate vector add both in sequential on CPU and in parallel on device.
//************************************
int main(int argc, char* argv[]) {
   
   // vectorOperation();
    //matrixProductVector();
    matrixInversion();

    std::cout << "Vector add successfully completed on device.\n";
    return 0;
}