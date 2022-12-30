//==============================================================
// Copyright © 2021 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#include "FemImpLibrary.h"

#include <iostream>
#include <iomanip>

#include "mkl.h"
#include "mkl_pblas.h"
#include "mkl_spblas.h"
#include <algorithm>
#include <functional>
#include <iterator>

#include <vector>
#include <string>
//using namespace sycl;

namespace
{
    template<class T>
    void sample(std::string& element, std::vector<T>& s, std::vector<T>& wt, int dim)
    {
        int nip = s.size() / dim;

        T root3 = T(1.0) / std::sqrt(T(3.0));
        T r15 = T(0.2) * std::sqrt(15.0);

        T w[3] = { T(5.0) / T(9.0), T(8.0) / T(9.0), T(5.0) / T(9.0) };
        T v[9] = { T(5.0) / T(9.0) * w[0],
                 T(5.0) / T(9.0) * w[1],
                 T(5.0) / T(9.0) * w[2],

                 T(8.0) / T(9.0) * w[0],
                 T(8.0) / T(9.0) * w[1],
                 T(8.0) / T(9.0) * w[2],

                 T(5.0) / T(9.0) * w[0],
                 T(5.0) / T(9.0) * w[1],
                 T(5.0) / T(9.0) * w[2] };

        if (element == "quadrilateral")
        {
            switch (nip)
            {
            case 9:
            {
                int i = 0;
                
                for (int j = 0; j < 7; j = j + 3)
                {
                    s[i + j * dim] = -r15;
                }

                for (int j = 1; j < 8; j = j + 3)
                {
                    s[i + j * dim] = 0.0;
                }

                for (int j = 2; j < nip; j = j + 3)
                {
                    s[i + j * dim] = r15;
                }
                                
                i = 1;

                for (int j = 0; j < 3; ++j)
                {
                    s[i + j * dim] = r15;
                }

                for (int j = 3; j < 6; ++j)
                {
                    s[i + j * dim] = 0.0;
                }

                for (int j = 6; j < nip; ++j)
                {
                    s[i + j * dim] = -r15;
                }

                std::copy(std::begin(v), std::end(v), wt.begin());
                break;
            }
            default:
                break;
            }
        }
        else if (element == "tetrahedron")
        {
            switch (nip)
            {
            case 1:
            {
                s[0] = 0.25; //s(1, 1) = 0.25_iwp
                s[1] = 0.25; //s(1, 2) = 0.25_iwp
                s[2] = 0.25; //s(1, 3) = 0.25_iwp
                wt[0] = 1.0 / 6.0;
            }
            default:
                break;
            }
        }
    }

    template<class T>
    void shape_der(std::vector<T>& der, const std::vector<T>& points, const int i, const int dim, const int nod)
    {
        T eta, xi, zeta, xi0, eta0, zeta0, etam, etap, xim, xip, c1, c2, c3;
        T t1, t2, t3, t4, t5, t6, t7, t8, t9, x2p1, x2m1, e2p1, e2m1, zetam, zetap;

        T   zero = 0.0, pt125 = 0.125, pt25 = 0.25, pt5 = 0.5,
            pt75 = 0.75, one = 1.0, two = 2.0, d3 = 3.0, d4 = 4.0, d5 = 5.0,
            d6 = 6.0, d8 = 8.0, d9 = 9.0, d10 = 10.0, d11 = 11.0,
            d12 = 12.0, d16 = 16.0, d18 = 18.0, d27 = 27.0, d32 = 32.0,
            d36 = 36.0, d54 = 54.0, d64 = 64.0, d128 = 128.0;

        int xii(20), etai(20), zetai(20), l;

        int nip = points.size() / dim;

        switch (dim)
        {
        case 2:
        {
            xi = points[i * dim];
            eta = points[1 + i * dim];
            c1 = xi;
            c2 = eta;
            c3 = one - c1 - c2;
            etam = pt25 * (one - eta);
            etap = pt25 * (one + eta);
            xim = pt25 * (one - xi);
            xip = pt25 * (one + xi);
            x2p1 = two * xi + one;
            x2m1 = two * xi - one;
            e2p1 = two * eta + one;
            e2m1 = two * eta - one;

            switch (nod)
            {
            case 8:
            {
                der[(1 - 1)] = etam * (two * xi + eta);
                der[(2 - 1)] = -d8 * etam * etap;
                der[(3 - 1)] = etap * (two * xi - eta);
                der[(4 - 1)] = -d4 * etap * xi;
                der[(5 - 1)] = etap * (two * xi + eta);
                der[(6 - 1)] = d8 * etap * etam;
                der[(7 - 1)] = etam * (two * xi - eta);
                der[(8 - 1)] = -d4 * etam * xi;
                der[(1 - 1) + nod] = xim * (xi + two * eta);
                der[(2 - 1) + nod] = -d4 * xim * eta;
                der[(3 - 1) + nod] = xim * (two * eta - xi);
                der[(4 - 1) + nod] = d8 * xim * xip;
                der[(5 - 1) + nod] = xip * (xi + two * eta);
                der[(6 - 1) + nod] = -d4 * xip * eta;
                der[(7 - 1) + nod] = xip * (two * eta - xi);
                der[(8 - 1) + nod] = -d8 * xim * xip;

                break;
            }

            default:
                break;

            }
            break;
        }
        default:
            break;
        }
    }

    template<class T>
    void shape_fun(std::vector<T>& fun, const std::vector<T>& points, const int i, const int dim, const int nod)
    {
        T eta, xi, etam, etap, xim, xip, zetam, zetap, c1, c2, c3;
        T t1, t2, t3, t4, t5, t6, t7, t8, t9;
        T zeta, xi0, eta0, zeta0;
        int xii(20), etai(20), zetai(20), l, ndim;
        T pt125 = 0.125, pt25 = 0.25, pt5 = 0.5, pt75 = 0.75,
            one = 1.0, two = 2.0, d3 = 3.0, d4 = 4.0, d8 = 8.0, d9 = 9.0,
            d16 = 16.0, d27 = 27.0, d32 = 32.0, d64 = 64.0, d128 = 128.0;
        
        switch (dim)
        {
        case 2:
        {
            c1 = points[i * dim];
            c2 = points[1 + i * dim];
            c3 = one - c1 - c2;
            xi = points[i * dim];
            eta = points[1 + i * dim];
            etam = pt25 * (one - eta);
            etap = pt25 * (one + eta);
            xim = pt25 * (one - xi);
            xip = pt25 * (one + xi);

            switch (nod)
            {
            case 8:
            {
                T funTEm[] = { d4 * etam * xim * (-xi - eta - one), d32 * etam * xim * etap,
                     d4 * etap * xim * (-xi + eta - one), d32 * xim * xip * etap,
                     d4 * etap * xip * (xi + eta - one), d32 * etap * xip * etam,
                     d4 * xip * etam * (xi - eta - one), d32 * xim * xip * etam };

                std::copy(std::begin(funTEm), std::end(funTEm), fun.begin());
            }
            default:
                break;
            }
            break;
        }
        default:
            break;
        }

    }

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

    template<class T>
    void deemat(T* dee, int nbRowsDee, T e, T v)
    {
        int ih = nbRowsDee;

        T v1, v2, c, vv, zero = 0.0, pt5 = 0.5, one = 1.0, two = 2.0;
      
        v1 = one - v;
        c = e / ((one + v) * (one - two * v));

        switch (ih)
        {
        case 3:
        {
            dee[0 + (0) *ih] = v1 * c;
            dee[1 + (1) *ih] = v1 * c;
            dee[2 + (2) *ih] = pt5 * c * (one - two * v);
            
            dee[1 + (0) * ih] = v * c;
            dee[0 + (1) * ih] = v * c;

            break;
        }

        case 6:
        {
            T v2 = v / (one - v);
            T vv = (one - two * v) / (one - v) * pt5;

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

            for (int i = 0; i < nbRowsDee * nbRowsDee; ++i)
            {
                dee[i] *= e / (two * (one + v) * vv);
            }

            break;
        }

        default:
            break;

        }
    }

    template<class T>
    void fsparv(T* kv, const T* km, int* g, const int* kdiag, const int idof)
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

    // 
    // This subroutine performs Cholesky factorisation on a symmetric
    // skyline global matrix.
    //
    template<class T>
    void sparin(T* kv, int* kdiag, int neq)
    {
        int n = neq;
        kv[0] = std::sqrt(kv[0]);

        for (int i = 2; i <= n; ++i)
        {
            T x = T(0.0);
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

    // 
    // This subroutine performs Cholesky forwardand back - substitution
    // on a symmetric skyline global matrix.
    //
    template<class T>
    void spabac(const T* kv, T* loads, const int* kdiag, int neq)
    {
        int n = neq;
        loads[1] = loads[1] / kv[0];

        for (int i = 2; i <= n; ++i)
        {
            int ki = kdiag[i - 1] - i;
            int l = kdiag[i - 1 - 1] - ki + 1;
            T x = loads[i];

            if (l != i)
            {
                int m = i - 1;
                for (int j = l; j <= m; ++j)
                {
                    x = x - kv[ki + j - 1] * loads[j];
                }
            }

            loads[i] = x / kv[ki + i - 1];

        }

        for (int it = 2; it <= n; ++it)
        {
            int i = n + 2 - it;
            int ki = kdiag[i - 1] - i;
            T x = loads[i] / kv[ki + i - 1];
            loads[i] = x;
            int l = kdiag[i - 1 - 1] - ki + 1;

            if (l != i)
            {
                int m = i - 1;
                for (int k = l; k <= m; ++k)
                {
                    loads[k] = loads[k] - x * kv[ki + k - 1];
                }

            }
        }

        loads[1] = loads[1] / kv[0];
    }
}

/// INTEL ///////////////////////////////
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

void copyVec(int n, const float* x, float* y)
{
    cblas_scopy(n, x, 1, y, 1);
}

void copyVec(int n, const double* x, double* y)
{
    cblas_dcopy(n, x, 1, y, 1);
}

void scalVecProduct(int n, const float a, float* x)
{
    cblas_sscal(n, a, x, 1);
}

void scalVecProduct(int n, const double a, double* x)
{
    cblas_dscal(n, a, x, 1);
}

float dotProduct(MKL_INT n, float* x, float* y)
{
    return cblas_sdot(n, x, 1, y, 1);
}

double dotProduct(MKL_INT n, double* x, double* y)
{
    return cblas_ddot(n, x, 1, y, 1);
}

void addVectors(const float* X, const float a, const float* Y, const float b, float* R, const float n)
{
    //x = a * x
    cblas_scopy(n, Y, 1, R, 1);
    cblas_sscal(n, b, R, 1);
    //y: = a * x + y
    cblas_saxpy(n, a, X, 1, R, 1);
}

void addVectors(const double* X, const double a, const double* Y, double b, double* R, const double n)
{
    //x = a * x
    cblas_dcopy(n, Y, 1, R, 1);
    cblas_dscal(n, b, R, 1);
    //y: = a * x + y
    cblas_daxpy(n, a, X, 1, R, 1);
}

void addVectors(const float a, const float* X, float* Y, const float n)
{
    //y: = a * x + y
    cblas_saxpy(n, a, X, 1, Y, 1);
}

void addVectors(const double a, const double* X, double* Y, const double n)
{
    //y: = a * x + y
    cblas_daxpy(n, a, X, 1, Y, 1);
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

void matricesAdd(const float* A, const float a, const float* B, const float b, float* C, const int m, const int n)
{
    mkl_somatadd('r', 'n', 'n', m, n, a, A, n, b, B, n, C, n);
}

void matricesAdd(const double* A, const double a, const double* B, const double b, double* C, const int m, const int n)
{
    mkl_domatadd('r', 'n', 'n', m, n, a, A, n, b, B, n, C, n);
}

// y = alpha * A * x + beta * y
void matmulVec(const float alpha, const float* A, const float* x, const float beta, float* y, const int m, const int n)
{
    cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, n, x, 1, beta, y, 1);
}

// y = alpha * A * x + beta * y
void matmulVec(const double alpha, const double* A, const double* x, const double beta, double* y, const int m, const int n)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, n, x, 1, beta, y, 1);
}

// C = alpha A * B + beta * C
void matmul(const float* A, const float* B, float* C, const int m, const int k, const int n)
{
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);
}

// C = alpha A * B + beta * C
void matmul(const double* A, const double* B, double* C, const int m, const int k, const int n)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);
}

void matmulTransA(const float* A, const float* B, float* C, const int m, const int k, const int n)
{
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, B, n, 0, C, n);
}

void matmulTransA(const double* A, const double* B, double* C, const int m, const int k, const int n)
{
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, B, n, 0, C, n);
}

// mkl_?getrfnp could be an alternative withouh using pivot
float invert(float* A, int m)
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

double invert(double* A, int m)
{
    lapack_int* ipiv = new lapack_int[m];

    auto info = LAPACKE_dgetrf(CblasRowMajor, m, m, A, m, ipiv);

    double determinant = 1.0;

    for (int i = 0; i < m; ++i)
    {
        determinant = determinant * A[i + i * m];
        determinant = ipiv[i] - 1 != i ? -determinant : determinant;
    }

    LAPACKE_dgetri(CblasRowMajor, m, A, m, ipiv);

    return determinant;
}

//////////////////////////////////////////////////////////////////////////////////

template<class T>
void beemat(T* bee, int nbRowsBee, int nbColumsBee, T* deriv, int nbColumsDerivs)
{
    int k, l, n, ih(nbRowsBee), nod(nbColumsDerivs);
    T x, y, z;

    switch (ih)
    {
    case 3:
    case 4:
    {
        for (int m = 1; m <= nod; ++m)
        {
            k = 2 * m;
            l = k - 1;
            x = deriv[m - 1 + 0 * nbColumsDerivs]; // deriv(1, m);
            y = deriv[m - 1 + 1 * nbColumsDerivs];
            bee[l - 1 + 0 * nbColumsBee] = x;
            bee[k - 1 + 2 * nbColumsBee] = x;
            bee[k - 1 + 1 * nbColumsBee] = y;
            bee[l - 1 + 2 * nbColumsBee] = y;
        }

        break;
    }
    case 6:
    {
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

        break;
    }
    default:
        break;
    }
    
}

namespace {
    
    template<class T>
    void ecmat(T* ecm, T* nt, T* tn, const T* fun, const int ndof, const int nodof)
    {
        int nod = ndof / nodof;
        
        for (int i = 1; i <= nod; ++i)
        {
            for (int j = 1; j <= nodof; ++j)
            {
                nt[((i - 1) * nodof + j - 1) * nodof + j - 1] = fun[i - 1];
                tn[(j - 1) * ndof + (i - 1) * nodof + j - 1] = fun[i - 1];
            }
        }

        matmul(nt, tn, ecm, ndof, nodof, ndof);
    }

    template<class T>
    void linmul_sky(const T* kv, const T* disps, T* loads, const int* kdiag, int neq)
    {
        //
        //This subroutine forms the product of symmetric matrix stored as
        //a skylineand a vector.
        //

        int n = neq;
        int low = 0;

        for (int i = 1; i <= n; ++i)
        {
            T x = T(0.0);
            int lup = kdiag[i - 1];
            if (i == 1)low = lup;
            if (i != 1)low = kdiag[i - 1 - 1] + 1;

            for (int j = low; j <= lup; ++j)
            {
                x = x + kv[j - 1] * disps[i + j - lup];
            }

            loads[i] = x;
            if (i != 1)
            {
                lup = lup - 1;
                for (int j = low; j <= lup; ++j)
                {
                    int k = i + j - lup - 1;
                    loads[k] = loads[k] + kv[j - 1] * disps[i];
                }
            }
        }
    }

    template<class T>
    void checon(const std::vector<T>& loads, std::vector<T>& oldIds, const T tol, bool& converged)
    {
        converged = true;

        T maxDiff = T(0.0);
        T maxLoads = T(0.0);

        for (int i = 0; i < loads.size(); ++i)
        {
            maxDiff = std::max(maxDiff, std::abs(loads[i] - oldIds[i]));
            maxLoads = std::max(maxLoads, std::abs(loads[i]));
        }

        converged = ((maxDiff / maxLoads) <= tol);

        std::copy(std::begin(loads), std::end(loads), oldIds.begin());
    }

    template<class T>
    float load(T t)
    {
        return std::sin(T(1.0) * t);
    }

    template<class T>
    void buildM(const std::vector<T>& m, std::vector<T>& M, int* num, int nod, int nodof, int nn)
    {
        //for (int i = 0; i < nod; ++i)
        //{
        //    for (int j = 0; j < nod; ++j)
        //    {
        //        if (num[i] <= num[j])
        //        {
        //            int K_Id_nn = num[j] * nodof + (num[i] * nodof * nn * nodof);
        //            int K_Id_nn_t = num[i] * nodof + (num[j] * nodof * nn * nodof);
        //
        //            int k_Id = j * nodof + i * nodof * nod * nodof;
        //            int k_Id_t = i * nodof + j * nodof * nod * nodof;
        //
        //            if (k_Id_t < k_Id)
        //                std::swap(k_Id_t, k_Id);
        //                                 
        //            for (int k = 0; k < nodof; ++k)
        //            {
        //                int KId = K_Id_nn + k;
        //                int KId_t = K_Id_nn_t + (k * nn * nodof);
        //                                        
        //                M[KId] += m[k_Id + k];
        //
        //                if (KId != KId_t)
        //                {
        //                    M[KId_t] = M[KId];
        //                }
        //            }
        //        }
        //    }
        //}

        for (int i = 0; i < nod * nodof; ++i)
        {
            for (int j = 0; j < nod * nodof; ++j)
            {
                int k_Id = j + i * nodof * nod;
                int K_Id = (num[j / nodof] * nodof + j % nodof) + ((num[i / nodof] * nodof + i % nodof) * nn * nodof);

                M[K_Id] += m[k_Id];
            }
        }
        //for (int i = 0; i < nn * nodof; i = i + 1)//nodof)
        //{
        //    for (int j = 0; j < nn * nodof; j = j + 1)//nodof)
        //    {
        //        int id = j + i * nn * nodof;
        //
        //        if (M[id] > 0)
        //        {
        //            std::cout << std::fixed;
        //            std::cout << std::setprecision(2);
        //            std::cout << "(" << M[id] << ")" << "  ";
        //        }
        //        else
        //        {
        //            std::cout << std::fixed;
        //            std::cout << std::setprecision(2);
        //            std::cout << M[id] << "  ";
        //        }
        //    }
        //    std::cout << std::endl;
        //    std::cout << std::endl;
        //}
        //
        //std::cout << "*******************" << std::endl;
        //std::cout << "*******************" << std::endl;
        //std::cout << "*******************" << std::endl;
    }

    template<class T>
    void buildSystem()
    {

    }
}

template<class T>
class Fem_Algoritm
{
private:
    //ndim = number of dimensions
    int ndim;

    //nodof = number of freedoms per node (x, y, z)
    int nodof;

    //nn = total number of nodes in the problem
    int nn;

    //nels = number of elements
    int nels;

    //neq = number of degree of freedom in the mesh
    int neq;

    //nod = number of node per element
    int nod;

    //nip = number of intregation points per element
    int nip;

    //nst = number of stress / strain terms
    int nst;

    //nf = nodal freedom array(nodof rows and nn colums)
    std::vector<int> nf;

    std::vector<T> points;
    std::vector<T> weights;

    T time = T(0);
    
    std::vector<T> fun;

    std::vector<T> der;

    std::vector<T> diag_precon;

    /////////////////////////
    std::vector<T> Km;
    std::vector<T> Mm;

    /////////////////////////

    std::vector<T> storkm;
    std::vector<T> stormm;

    //cg
    std::vector<T> A;
    //

    //std::vector<T> u;
    //std::vector<T> p;
    //std::vector<T> xnew;
    
    std::vector<int> node;
    std::vector<T> val;


    // skyline profile
    //std::vector<int> kdiag;

    //kv = global stiffness matrix
    //std::vector<T> kv;
    std::vector<int> g_g;
    //global consisten mass
    //std::vector<T> mv;

    //left hand side matrix (stored as a skyline) 
    std::vector<T> f1;

    //gravlo = global gravity loading vector 
    std::vector<T> gravlo;
    T gravity = T(980);
    std::vector<T> gravityDir;

    std::vector<T> x0;
    std::vector<T> d1x0;
    std::vector<T> x1;
    std::vector<T> d2x0;
    std::vector<T> d1x1;
    std::vector<T> d2x1;


    std::vector<float> x02;
    std::vector<float> d1x02;
    std::vector<float> x12;
    std::vector<float> d2x02;
    std::vector<float> d1x12;
    std::vector<float> d2x12;

    ///////////////////////////////////
    std::vector<T> prop;
    T fm = T(0.0001);
    T fk = T(0.0);
    ///////////////////////////////////
    // eld = element displacement
    std::vector<T> eld;
    std::vector<T> loads;

    int itt = 1;

    std::string   element;

public:
    
    void setGravityAcceleration(T gravAcc)
    {
        gravity = gravAcc;
    }

    void setGravityDirection(T* gravDir)
    {
        gravityDir[0] = gravDir[0];
        gravityDir[1] = gravDir[1];
        gravityDir[2] = gravDir[2];
    }

    void setDamping(const T fk, const T fm)
    {
        this->fm = fm;
        this->fk = fk;
    }

    void setMaterialParams(const T e, const T v, const T gamma)
    {
        if (!prop.empty())
        {
            prop[0] = e;
            prop[1] = v;
            prop[2] = gamma;
        }
    }

    void loadedNodes(int* nodes, int loaded_nodes, T* vals)
    {
        node.resize(loaded_nodes);
        val.resize(loaded_nodes * ndim);
        std::copy(nodes, nodes + loaded_nodes, node.begin());
        std::copy(vals, vals + loaded_nodes * ndim, val.begin());
    }

public:
    Fem_Algoritm(int ndim, int nodof, int nels, int nod, int nip, const char* element)
        : ndim(ndim), nodof(nodof), nod(nod), element(element), nip(nip), nn(0), nels(nels), neq(0)
    {
        fun.resize(nod, T(0.25));

        if (std::string("tetrahedron") == std::string(element))
        {
            T derTemp[4 * 3] =
            {
                    T(1.0), T(0.0), T(0.0), T(-1.0),
                    T(0.0), T(1.0), T(0.0), T(-1.0),
                    T(0.0), T(0.0), T(1.0), T(-1.0)
            };
            der.insert(der.end(), std::begin(derTemp), std::end(derTemp));

            nst = 6;
        }
        else
        {
            der.resize(ndim * nod, T(0.0));
        }

        if (std::string("quadrilateral") == std::string(element))
        {
            nst = 3;
        }

        points.resize(nip * ndim, T(0.25));
        weights.resize(nip, T(1.0) / T(6.0));

        //nprops = number of material properties
        const int nprops = 3;

        //np_types = number of diffent property types
        const int np_types = 1;

        //prop = material property(e, v, gamma)
        prop.resize(nprops * np_types, 0.0f);

        prop[0] = T(100000000.0);
        prop[1] = T(0.3);
        prop[2] = T(333);

        //
        gravityDir.resize(ndim, 0.0);
        gravityDir[ndim - 1] = T(1.0);
    }

    void init(T* g_coord, int* g_num, int* in_nf, int in_nn)
    {
        nn = in_nn;
               
        nf.resize(nodof * nn, 0);
        
        std::copy(in_nf, in_nf + nodof * nn, std::begin(nf));
               
        neq = formnf(&nf[0], nodof, nn);
               
        //ndof = number of degree of freedom per element
        int ndof = nod * nodof;
                 
        g_g.resize(ndof * nels, 0);
                
        /////////////////
        std::vector<float> g_coord_T(nodof * nn, 0.0f);

        for (int k = 0; k < nn; ++k)
        {
            for (int l = 0; l < nodof; ++l)
            {
                g_coord_T[k + l * nn] = g_coord[l + k * nodof];
            }
        }
        /////////////////
 
        //---------------------- loop the elements to find global arrays sizes---- -
                
        for (int i = 0; i < nels; ++i)
        {
            std::vector<int> g(ndof, 0);
            int* num = &g_num[nod * i];
        
            num_to_g(num, &nf[0], &g[0], nod, nodof, nn);
        
            for (int j = 0; j < ndof; ++j)
            {
                g_g[i + j * nels] = g[j];
            }
        }
      
        sample(element, points, weights, ndim);
        gravlo.resize(neq + 1, T(0.0));
        diag_precon.resize(neq + 1, T(0.0));

        //----element stiffnessand mass integration, storage and preconditioner-- -

        const T e = prop[0];
        const T v = prop[1];
                        
        eld.resize(ndof, 0.0);
        
        ///
        Km.resize(nodof * nn * nodof * nn, T(0.0f));
        Mm.resize(nodof * nn * nodof * nn, T(0.0f));
        ///

        std::vector<T> verify;

        for (int i = 0; i < nels; ++i)
        {
            std::vector<T> dee(nst * nst, T(0.0));
            deemat(&dee[0], nst, e, v);
        
            int* num = &g_num[nod * i];
        
            std::vector<T> coord(nod * ndim, T(0.0));
        
            for (int j = 0; j < nod; ++j)
            {
                for (int k = 0; k < ndim; ++k)
                {
                    coord[k + j * ndim] = g_coord_T[num[j] + k * nn];
                }
            }
                
            std::vector<T> km(ndof * ndof, T(0.0));
            std::vector<T> mm(ndof * ndof, T(0.0));
            ////// for each point of integration, we have just one for 3d tets
            
            for (int j = 0; j < nip; ++j)
            {
                std::vector<T> jac(ndim * ndim, T(0.0));
                std::vector<T> deriv(ndim * nod, T(0.0));
                std::vector<T> bee(nst * ndof, T(0.0));

                shape_der(der, points, j, ndim, nod);
                shape_fun(fun, points, j, ndim, nod);

                // calculate jac
                matmul(&der[0], &coord[0], &jac[0], ndim, nod, ndim);
                T det = invert(&jac[0], ndim);

                // calculate the derivative in x, y, z
                matmul(&jac[0], &der[0], &deriv[0], ndim, ndim, nod);

                std::fill(std::begin(bee), std::end(bee), T(0.0));
                beemat(&bee[0], nst, ndof, &deriv[0], nod);

                std::vector<T> temporal(nst * ndof, T(0.0));
                std::vector<T> kmTemp(ndof * ndof, T(0.0));

                matmulTransA(&bee[0], &dee[0], &temporal[0], ndof, nst, nst);
                matmul(&temporal[0], &bee[0], &kmTemp[0], ndof, nst, ndof);

                std::for_each(std::begin(kmTemp), std::end(kmTemp), [&](T& v) { v *= det * weights[j]; });

                matricesAdd(&km[0], T(1.0), &kmTemp[0], T(1.0), &km[0], ndof, ndof);

                std::vector<T> ecm(ndof * ndof, T(0.0));
                std::vector<T> nt(ndof * nodof, T(0.0));
                std::vector<T> tn(nodof * ndof, T(0.0));

                ecmat(&ecm[0], &nt[0], &tn[0], &fun[0], ndof, nodof);

                int etype = 1;
                T bEcm = det * weights[j] * prop[2 * etype];
                matricesAdd(&mm[0], T(1.0), &ecm[0], bEcm, &mm[0], ndof, ndof);

                //int kk = 0;
                //for (int k = nodof; k <= ndof; k = k + nodof) // TODO: generalize to a general displacement.
                //{
                //    eld[k - 1] = eld[k - 1] + fun[kk] * det * weights[j];
                //    eld[k - 1 + 1] = eld[k - 1 + 1] + fun[kk] * det * weights[j];
                //    eld[k - 1 + 2] = eld[k - 1 + 2] + fun[kk] * det * weights[j];
                //
                //    kk++;
                //}
            }
            
            buildM(mm, Mm, num, nod, nodof, nn);
            buildM(km, Km, num, nod, nodof, nn);

            storkm.insert(storkm.end(), km.begin(), km.end());
            stormm.insert(stormm.end(), mm.begin(), mm.end());
        }
              
        A.resize(storkm.size(), T(0.0));

        //---------------------- - initial conditions and factorise equations--------
        x0.resize(neq + 1, T(0.0));
        x1.resize(neq + 1, T(0.0));

        for (int i = 0; i < nn; ++i)
        {
            x0[nf[i]] = (nf[i] > 0) ? g_coord[i * 3] : T(0.0);
            x0[nf[i + nn]] = (nf[i + nn] > 0) ? g_coord[i * 3 + 1] : T(0.0);
            x0[nf[i + 2 * nn]] = (nf[i + 2 * nn] > 0) ? g_coord[i * 3 + 2] : T(0.0);
        }

        std::copy(x0.begin(), x0.end(), x1.begin());
               

        d1x0.resize(neq + 1, T(0.0));
        d2x0.resize(neq + 1, T(0.0));
        d1x1.resize(neq + 1, T(0.0));
        d2x1.resize(neq + 1, T(0.0));

        ////////////////////////////////////////////////////////////

        x02.resize(neq, T(0.0));
        x12.resize(neq, T(0.0));

        for (int i = 0; i < nn; ++i)
        {
            if (nf[i] > 0)
                x02[nf[i] - 1] = g_coord[i * 3];

            
            if (nf[i + nn] > 0)
                x02[nf[i + nn] - 1] = g_coord[i * 3 + 1];

            if (nf[i + 2 * nn] > 0)
                x02[nf[i + 2 * nn] - 1] = g_coord[i * 3 + 2];
        }

        std::copy(x02.begin(), x02.end(), x12.begin());


        d1x02.resize(neq, T(0.0));
        d2x02.resize(neq, T(0.0));
        d1x12.resize(neq, T(0.0));
        d2x12.resize(neq, T(0.0));

        ////////////////////////////////////////////////////////////


        loads.resize(neq + 1, T(0.0));
    }
    
    void cg(std::vector<T>& xnew)
    {
        int ndof = nod * nodof;
        int cg_iters = 0;
        int cg_limit = 50;
        T cg_tol = T(0.00001);

        std::vector<T> d(neq + 1, T(0.0));

        for (int i = 0; i < neq + 1; ++i)
        {
            d[i] = diag_precon[i] * loads[i];
        }

        std::vector<T> p(d);
        std::vector<T> x(neq + 1, T(0.0));
        std::vector<T> u(loads);
        
        bool cg_converged(false);
        do
        {
            cg_iters = cg_iters + 1;
            std::fill(std::begin(u), std::end(u), T(0.0));

            for (int i = 0; i < nels; ++i)
            {
                std::vector<T> uTemp(ndof, 0);
                std::vector<T> pTemp(ndof, 0);

                std::vector<int> g(ndof, 0);
                for (int j = 0; j < ndof; ++j)
                {
                    int g_index = g_g[i + j * nels];

                    pTemp[j] = p[g_index];
                    g[j] = g_index;
                }

                const T* aij = &A[i * ndof * ndof];
                matmulVec(T(1.0), aij, &pTemp[0], T(1.0), &uTemp[0], ndof, ndof);

                for (int j = 0; j < ndof; ++j)
                {
                    int g_index = g_g[i + j * nels];

                    u[g_index] += uTemp[j];
                }
            }

            u[0] = T(0.0);

            T up = dotProduct(neq + 1, &loads[0], &d[0]);
            T alpha = up / dotProduct(neq + 1, &p[0], &u[0]);

            addVectors(&x[0], T(1.0), &p[0], alpha, &xnew[0], neq + 1); // xnew = x + alpha * p
            addVectors(-alpha, &u[0], &loads[0], neq + 1); // loads = loads - alpha * u

            for (int i = 0; i < neq + 1; ++i)
            {
                d[i] = diag_precon[i] * loads[i];
            }

            T beta = dotProduct(neq + 1, &loads[0], &d[0]) / up;
            addVectors(1 / beta, &d[0], &p[0], neq + 1);
            scalVecProduct(neq + 1, beta, &p[0]);

            checon(xnew, x, cg_tol, cg_converged);

            if (cg_converged)
            {
                break;
            }
        } while (cg_iters < cg_limit);
    }

    void cg2(sparse_matrix_t& csrA, matrix_descr& descrA, std::vector<float>& xnew, std::vector<float>&  diag_precon2, std::vector<float>& loads2)
    {
        int ndof = nod * nodof;
        int cg_iters = 0;
        int cg_limit = 50;
        float cg_tol = float(0.00001);

        std::vector<float> d(neq, float(0.0));

        for (int i = 0; i < neq; ++i)
        {
            d[i] = diag_precon2[i] * loads2[i];
        }

        std::vector<float> p(d);
        std::vector<float> x(neq, float(0.0));
        std::vector<float> u(loads2);

        bool cg_converged(false);
        do
        {
            cg_iters = cg_iters + 1;
            std::fill(std::begin(u), std::end(u), float(0.0));

            //matmulVec(T(1.0), aij, &pTemp[0], T(1.0), &uTemp[0], ndof, ndof);
            
            sparse_status_t status = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, csrA, descrA, &p[0], 1.0f, &u[0]);

            //for (int i = 0; i < nels; ++i)
            //{
            //    std::vector<T> uTemp(ndof, 0);
            //    std::vector<T> pTemp(ndof, 0);
            //
            //    std::vector<int> g(ndof, 0);
            //    for (int j = 0; j < ndof; ++j)
            //    {
            //        int g_index = g_g[i + j * nels];
            //
            //        pTemp[j] = p[g_index];
            //        g[j] = g_index;
            //    }
            //
            //    const T* aij = &A[i * ndof * ndof];
            //    matmulVec(T(1.0), aij, &pTemp[0], T(1.0), &uTemp[0], ndof, ndof);
            //
            //    for (int j = 0; j < ndof; ++j)
            //    {
            //        int g_index = g_g[i + j * nels];
            //
            //        u[g_index] += uTemp[j];
            //    }
            //}

            //u[0] = T(0.0);

            float up = dotProduct(neq, &loads2[0], &d[0]);
            float alpha = up / dotProduct(neq, &p[0], &u[0]);

            addVectors(&x[0], T(1.0), &p[0], alpha, &xnew[0], neq); // xnew = x + alpha * p
            addVectors(-alpha, &u[0], &loads2[0], neq); // loads = loads - alpha * u

            for (int i = 0; i < neq; ++i)
            {
                d[i] = diag_precon2[i] * loads2[i];
            }

            float beta = dotProduct(neq, &loads2[0], &d[0]) / up;
            addVectors(1 / beta, &d[0], &p[0], neq);
            scalVecProduct(neq, beta, &p[0]);

            checon<float>(xnew, x, cg_tol, cg_converged);

            if (cg_converged)
            {
                break;
            }
        } while (cg_iters < cg_limit);
    }

    void update(T dtim, T* verticesBuffer)
    {
        //std::cout << "update 1" << std::endl;
        //// for debug purposes //
        time = time + dtim;

        int ndof = nod * nodof;
        int etype = 1;

        std::fill(gravlo.begin(), gravlo.end(), T(0.0));
        std::fill(A.begin(), A.end(), T(0.0));

        for (int i = 0; i < nels; ++i)
        {
            std::vector<int> g(ndof, 0);
            for (int j = 0; j < ndof; ++j)
            {
                g[j] = g_g[i + j * nels];
            }

            for (int k = 0; k < ndof; ++k)
            {
                gravlo[g[k]] = gravlo[g[k]] + prop[2 * etype] * gravityDir[k % ndim] * gravity;  //eld[k] *
            }
        }
        
        
        //////////////////////////
        //
        T c1 = fm;
        T c2 = fk;
        T c3 = T(1.0) + fm * dtim;
        T c4 = dtim + fk;
        
        //////////////////////////////////////////////////////////////

        //build system A and u
        std::fill(std::begin(loads), std::end(loads), T(0.0));

        for (int i = 0; i < nels; ++i)
        {
            std::vector<T> dxTemp(ndof, 0);
            std::vector<T> d1x0Temp(ndof, 0);

            for (int j = 0; j < ndof; ++j)
            {
                int g_index = g_g[i + j * nels];

                dxTemp[j] = x1[g_index] - x0[g_index];
                d1x0Temp[j] = d1x0[g_index];
            }

            const T* km = &storkm[i * ndof * ndof];
            const T* mm = &stormm[i * ndof * ndof];
            T* aij = &A[i * ndof * ndof];
                        
            matricesAdd(km, c4 * dtim, mm, c1, aij, ndof, ndof);

            /////////
            //u(g) = u(g) + MATMUL(km, dxTemp(g)) + MATMUL(km * c4 + mm * c1, d1x0(g)) // solving by acceleration
            std::vector<T> uTemp(ndof, 0);
            matmulVec(c4, km, &d1x0Temp[0], T(1.0), &uTemp[0], ndof, ndof);
            matmulVec(c1, mm, &d1x0Temp[0], T(1.0), &uTemp[0], ndof, ndof);

            matmulVec(T(1.0), km, &dxTemp[0], T(1.0), &uTemp[0], ndof, ndof);
            /////////
  
            for (int k = 0; k < ndof; ++k)
            {
                int g_index = g_g[i + k * nels];
                loads[g_index] += (T(-1.0) * uTemp[k]);
            }
        }
        
        loads[0] = T(0.0);

        //////////////////////////////////////////////////////////////

        //int loaded_nodes = node.size();
        //            
        //T temporal = dtim * load(time);
        //for (int j = 1; j <= loaded_nodes; ++j)
        //{
        //    for (int k = 0; k < ndim; ++k)
        //    {
        //        loads[nf[nn * k + node[j - 1] - 1]] +=
        //            val[(j - 1) * loaded_nodes + k] * temporal;
        //    }
        //}
        // 
        
        //////////////////////////////////////////////////////////////

        // build precondition
        for (int i = 0; i < nels; ++i)
        {
            std::vector<int> g(ndof, 0);
            for (int j = 0; j < ndof; ++j)
            {
                g[j] = g_g[i + j * nels];
            }

            for (int k = 0; k < ndof; ++k)
            {
                diag_precon[g[k]] = A[k + k * ndof];
            }
        }

        for (int k = 1; k < diag_precon.size(); ++k)
        {
            diag_precon[k] = T(1.0) / diag_precon[k];
        }
        diag_precon[0] = T(0.0);

        //////////////////////////////////////////////////////////////
                
        addVectors(T(1.0), &gravlo[0], &loads[0], loads.size());
                
        std::vector<T> xnew(neq + 1, T(0.0));

        //---------------------- - pcg equation solution----------------------------
        cg(xnew);
        //--------------------------------------------------------------------------
        
        //d1x1 = d1x0 + xnew * dtim;
        addVectors(&d1x0[0], T(1.0), &xnew[0], dtim, &d1x1[0], neq + 1);
                        
        //x1 = x1 + d1x1 * dtim;
        addVectors(dtim, &d1x1[0], &x1[0], neq + 1);
      
        //addVectors(&x1[0], T(1.0), &d1x1[0], dtim, &x1[0], neq + 1);
        
        copyVec(neq + 1, &d1x1[0], &d1x0[0]);
     


        ///////////////////
        //if (!node.empty())
        //{
        //    int npri = 1;
        //    static itt = 1;
        //    int nres = node[0];
        //
        //    if (itt / npri * npri == itt)
        //    {
        //        std::cout << "time " << time << "      load  " << load(time) << " node " << nres << "     x " << x1[nf[nres - 1]] << "     y " << x1[nf[nres + nn - 1]] << " cg it " << cg_iters << std::endl;// "     z obj   " << x0[nf[nres + 2 * nn - 1]] << std::endl;
        //    }
        //
        //    itt = itt + 1;
        //}
        ////////////////////

        if (verticesBuffer)
        {
            if (ndim == 3)
            {
                for (int i = 0; i < nn; ++i)
                {
                    verticesBuffer[i * 3] = (nf[i] > 0) ? x1[nf[i]]                       : verticesBuffer[i * 3];
                    verticesBuffer[i * 3 + 1] = (nf[i + nn] > 0) ? x1[nf[i + nn]]         : verticesBuffer[i * 3 + 1];
                    verticesBuffer[i * 3 + 2] = (nf[i + 2 * nn] > 0) ? x1[nf[i + 2 * nn]] : verticesBuffer[i * 3 + 2];

                    //if (i == 5)
                    //{
                    //    std::cout << x0[nf[i]] << " " << x0[nf[i + nn]] << " " << x0[nf[i + 2 * nn]] << std::endl;
                    //}
                }
            }
        }
    }

    void update2(T dtim, T* verticesBuffer)
    {
        //std::cout << "update 2" << std::endl;
        //// for debug purposes //
        time = time + dtim;

        int ndof = nod * nodof;
        int etype = 1;

        std::vector<float> gravlo2(neq, float(0.0));
        //std::fill(gravlo2.begin(), gravlo2.end(), T(0.0));
        
        
        for (int i = 0; i < nels; ++i)
        {
            std::vector<int> g(ndof, 0);
            for (int j = 0; j < ndof; ++j)
            {
                g[j] = g_g[i + j * nels];
            }

            for (int k = 0; k < ndof; ++k)
            {
                if (g[k] > 0)
                    gravlo2[g[k] - 1] = gravlo2[g[k] - 1] + prop[2 * etype] * gravityDir[k % ndim] * gravity;  //eld[k] *
            }
        }


        //////////////////////////
        //
        T c1 = fm;
        T c2 = fk;
        T c3 = T(1.0) + fm * dtim;
        T c4 = dtim + fk;

        //////////////////////////////////////////////////////////////

        //build system A and u
       /////////////////////////////////////////////////////////////


        std::vector<int> tempnf(nodof * nn, 0);
        
        std::vector<float> SPvalue_dt; //(neq, T(0.0));
        std::vector<float> SPvalue_Km; //(neq, T(0.0));
        std::vector<float> SPvalue; //(neq, T(0.0));
        std::vector<long long> SPcolumns;
        std::vector<long long> SProwIndex;// (neq + 1, 0);

        for (int k = 0; k < nn; ++k)
        {
            for (int l = 0; l < nodof; ++l)
            {
                tempnf[l + k * nodof] = nf[k + l * nn];
            }
        }

        //std::cout << "********************************" << std::endl;
        //std::cout << "********************************" << std::endl;
        //std::cout << "********************************" << std::endl;

        //for (int i = 0; i < nn * nodof; ++i)
        //{
        //    for (int j = 0; j < nn * nodof; ++j)
        //    {
        //        int id = j + i * nn * nodof;
        //
        //        float aij = Km[id] * c4 * dtim + Mm[id] * c1;
        //        //float aij = Km[id];
        //        //float aij = Mm[id];
        //
        //
        //        if (i <= j)
        //        {
        //            if (tempnf[j] > 0 && (tempnf[i] > 0))//std::abs(aij) > 0.0001)
        //            {
        //                std::cout << std::fixed;
        //                std::cout << std::setprecision(2);
        //                std::cout << "(" << aij << ")" << "  ";
        //            }
        //            else
        //            {
        //                std::cout << std::fixed;
        //                std::cout << std::setprecision(2);
        //                std::cout << aij << "  ";
        //            }
        //        }
        //    }
        //
        //    std::cout << std::endl;
        //    std::cout << std::endl;
        //}
                
        std::vector<float> dxTemp(neq, float(0.0));
        std::vector<float> d1x0Temp(neq, float(0.0));
        std::vector<float> diag_precon2;

        for (int i = 0; i < nn * nodof; ++i)
        {
            for (int j = 0; j < nn * nodof; ++j)
            {
                int id = j + i * nn * nodof;
                float aij_dt = Km[id] * c4 * dtim + Mm[id] * c1;
                float aij = Km[id] * c4 + Mm[id] * c1;

                if (tempnf[j] > 0 && (tempnf[i] > 0))
                {
                    dxTemp[tempnf[i] - 1] = x12[tempnf[i] - 1] - x02[tempnf[i] - 1];
                    d1x0Temp[tempnf[i] - 1] = d1x02[tempnf[i] - 1];

                    if (tempnf[i] <= tempnf[j])
                    {
                        if (tempnf[j] == tempnf[i])
                        {
                            SPvalue.push_back(aij);
                            SPvalue_dt.push_back(aij_dt);
                            SPvalue_Km.push_back(Km[id]);
                            SPcolumns.push_back(tempnf[j] - 1);
                            SProwIndex.push_back(SPvalue.size() - 1);
                            diag_precon2.push_back(aij);
                        }
                        else if (std::abs(aij_dt) > 0.00001)
                        {
                            SPvalue.push_back(aij);
                            SPvalue_dt.push_back(aij_dt);
                            SPvalue_Km.push_back(Km[id]);
                            SPcolumns.push_back(tempnf[j] - 1);
                        }
                    }
                }
            }
        }

        SProwIndex.push_back(SPvalue_dt.size());

        sparse_matrix_t csrA;
        sparse_matrix_t csrA_dt;
        sparse_matrix_t csrA_Km;

        sparse_status_t status;

        status = mkl_sparse_s_create_csr(&csrA,
            SPARSE_INDEX_BASE_ZERO,
            neq,  // number of rows
            neq,  // number of cols
            &SProwIndex[0],
            &SProwIndex[0] + 1,
            &SPcolumns[0],
            &SPvalue[0]);

        status = mkl_sparse_s_create_csr(&csrA_dt,
            SPARSE_INDEX_BASE_ZERO,
            neq,  // number of rows
            neq,  // number of cols
            &SProwIndex[0],
            &SProwIndex[0] + 1,
            &SPcolumns[0],
            &SPvalue_dt[0]);

        status = mkl_sparse_s_create_csr(&csrA_Km,
            SPARSE_INDEX_BASE_ZERO,
            neq,  // number of rows
            neq,  // number of cols
            &SProwIndex[0],
            &SProwIndex[0] + 1,
            &SPcolumns[0],
            &SPvalue_Km[0]);

        //////////////////////////////////////////////////////////////

        //u(g) = u(g) + MATMUL(km, dxTemp(g)) + MATMUL(km * c4 + mm * c1, d1x0(g)) // solving by acceleration

        std::vector<float> uTemp(neq, float(0.0));

        matrix_descr descrA;
        descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        descrA.mode = SPARSE_FILL_MODE_UPPER;
        descrA.diag = SPARSE_DIAG_NON_UNIT;
        status = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, csrA_Km, descrA, &dxTemp[0], 0.0f, &uTemp[0]);
        status = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, csrA, descrA, &d1x0Temp[0], 1.0f, &uTemp[0]);

        std::vector<float> loads2(neq, 0.0f);

        //        loads[g_index] += (T(-1.0) * uTemp[k]);

        //y: = a * x + y

        addVectors(float(-1.0), &uTemp[0], &loads2[0], neq);

        for (int k = 0; k < diag_precon2.size(); ++k)
        {
            diag_precon2[k] = float(1.0) / diag_precon2[k];
        }

        addVectors(float(1.0), &gravlo2[0], &loads2[0], loads2.size());
        
        
        std::vector<float> xnew(neq, float(0.0));
        //
        ////---------------------- - pcg equation solution----------------------------
        //cg(xnew);
        cg2(csrA_dt, descrA, xnew, diag_precon2, loads2);
        ////--------------------------------------------------------------------------
        
        //d1x1 = d1x0 + xnew * dtim;
        addVectors(&d1x02[0], float(1.0), &xnew[0], dtim, &d1x12[0], neq);
        
        //x1 = x1 + d1x1 * dtim;
        addVectors(dtim, &d1x12[0], &x12[0], neq);
        
        //addVectors(&x1[0], T(1.0), &d1x1[0], dtim, &x1[0], neq + 1);
        
        copyVec(neq, &d1x12[0], &d1x02[0]);

        if (verticesBuffer)
        {
            if (ndim == 3)
            {
                for (int i = 0; i < nn; ++i)
                {
                    if (nf[i] > 0)
                        verticesBuffer[i * 3] = x12[nf[i] - 1];

                    if (nf[i + nn] > 0)
                        verticesBuffer[i * 3 + 1] = x12[nf[i + nn] - 1];

                    if (nf[i + 2 * nn] > 0)
                        verticesBuffer[i * 3 + 2] = x12[nf[i + 2 * nn] - 1];
                                                
                    //if (i == 5)
                    //{
                    //    std::cout << x0[nf[i]] << " " << x0[nf[i + nn]] << " " << x0[nf[i + 2 * nn]] << std::endl;
                    //}
                }
            }
        }

        mkl_sparse_destroy(csrA);
        mkl_sparse_destroy(csrA_dt);
        mkl_sparse_destroy(csrA_Km);
    }
};

template<class T>
std::vector<Fem_Algoritm<T>*> FEM_Factory<T>::femAlg;

template<class T>
int FEM_Factory<T>::create(int ndim, int nodof, int nels, int nod, int nip, const char* element)
{
    femAlg.push_back(new Fem_Algoritm<T>(ndim, nodof, nels, nod, nip, element));
    return femAlg.size() - 1;
}

template<class T>
void FEM_Factory<T>::init(int id, T* g_coord, int* g_num, int* in_nf, int in_nn)
{
    if (id < femAlg.size())
    {
        femAlg[id]->init(g_coord, g_num, in_nf, in_nn);
    }
}

template<class T>
void FEM_Factory<T>::loadedNodes(int id, int* nodes, int loaded_nodes, T* vals)
{
    if (id < femAlg.size())
    {
        femAlg[id]->loadedNodes(nodes, loaded_nodes, vals);
    }
}

template<class T>
void FEM_Factory<T>::setGravityAcceleration(int id, T gravAcc)
{
    if (id < femAlg.size())
    {
        femAlg[id]->setGravityAcceleration(gravAcc);
    }
}

template<class T>
void FEM_Factory<T>::setGravityDirection(int id, T* gravDir)
{
    if (id < femAlg.size())
    {
        femAlg[id]->setGravityDirection(gravDir);
    }
}

template<class T>
void FEM_Factory<T>::setDamping(int id, const T fk, const T fm)
{
    if (id < femAlg.size())
    {
        femAlg[id]->setDamping(fk, fm);
    }
}

template<class T>
void FEM_Factory<T>::setMaterialParams(int id, const T e, const T v, const T gamma)
{
    if (id < femAlg.size())
    {
        femAlg[id]->setMaterialParams(e, v, gamma);
    }
}


template<class T>
void FEM_Factory<T>::update(int id, T dtim, T* verticesBuffer)
{
    if (id < femAlg.size())
    {
        //femAlg[id]->update2(dtim, verticesBuffer);
        femAlg[id]->update(dtim, verticesBuffer);
    }
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

float getValueFromExportSparse(int r, int c, MKL_INT* pointerB_C, MKL_INT* pointerE_C, MKL_INT* columns_C, float* values_C)
{
    //1.
    int accumulate = 0;

    for (int i = 0; i < r; ++i)
    {
        accumulate += pointerE_C[i] - pointerB_C[i];
    }

    //2.
    int first = accumulate;
    int second = accumulate + pointerE_C[r] - pointerB_C[r];
    auto itt = std::find(&columns_C[first], &columns_C[second], c);

    //3.
    auto pos = std::distance(&columns_C[0], itt);

    return values_C[pos];
}

FEMIMP_DLL_API void testSparse()
{
    const int m = 6;

    float matrix[m * m] = {
         100,    0.0,   0.0,  -3.5, -1.0,  0.8,
         0.0,    0.0,   0.0,   0.0,  0.0,  0.0, 
         0.0,    0.0,   4.0,  -2.1,  0.0,  6.1,
        -3.5,    0.0,  -2.1,  -3.8,  55.0, 0.0, 
        -1.0,    0.0,   0.0,   55.0, 7.0,  0.0,
         0.8,    0.0,   6.1,   0.0,  0.0,  -1.4
    };

    float vec[] = { 2.2, -3.3, 1, 4.5, 0.1, -0.8 };
    std::vector<float> result(m, 0.0f);
    matmulVec(1.0, matrix, vec, 0.0, &result[0], m, m);
    
    // Symetric Matrix CSR from m x m matrix
    // value[N]: upper triangle or lower triangle, storage from left to right
    // columns[N]: index column of the non-zero values, from left to right
    // rowIndex[m + 1]: First non-zero value (index value table) from left to right 
    //                 (the last value) is the number of the non-zero values.

    std::vector<float> value;
    std::vector<long long> columns;
    std::vector<long long> rowIndex;

    for (int i = 0; i < m; ++i)
    {
        for(int j = i; j < m; ++j)
        {
            int index = j + i * m;

            float v = matrix[index];

            if (i == j)
            {
                value.push_back(matrix[index]);
                columns.push_back(j);
                rowIndex.push_back(value.size() - 1);
            }
            else if (std::abs(v) > 0.000001)
            {
                value.push_back(matrix[index]);
                columns.push_back(j);
            }
        }
    }

    rowIndex.push_back(value.size());

    ////////////////
    sparse_matrix_t csrA;
    sparse_status_t status;

    status = mkl_sparse_s_create_csr(&csrA,
        SPARSE_INDEX_BASE_ZERO,
        m,  // number of rows
        m,  // number of cols
        &rowIndex[0],
        &rowIndex[0] + 1,
        &columns[0],
        &value[0]);
            
    std::vector<float> result2(m, 0.0f);

    matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    status = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, csrA, descrA, vec, 0.0, &result2[0]);

    for (int i = 0; i < m; ++i)
    {
        std::cout << " result with sgemv " << result[i] << " result with sparse_s_mv " << result2[i] << std::endl;
    }

    ////////
    // export
    sparse_index_base_t    indexing;
    struct matrix_descr    descr_type_gen;
    MKL_INT  rows, cols;
    MKL_INT* pointerB_C = NULL, * pointerE_C = NULL;
    MKL_INT *columns_C = NULL;
    float* values_C;

    status = mkl_sparse_s_export_csr( csrA,
                                      &indexing,
                                      &rows,
                                      &cols,
                                      &pointerB_C,
                                      &pointerE_C,
                                      &columns_C,
                                      &values_C);


    //
    int r = 3;
    int c = 4;

    float valueResq = getValueFromExportSparse(r, c, pointerB_C, pointerE_C, columns_C, values_C);
    std::cout << "Find val. Row  " << r << " Colum " << c << " value " << valueResq << std::endl;
    //

    ////////
    // set value
    float valToSet = 25.0f;
    status = mkl_sparse_s_set_value(csrA, r, c, valToSet);
    
    // verification//////////////
    /////////////////////////////

    status = mkl_sparse_s_export_csr(csrA,
        &indexing,
        &rows,
        &cols,
        &pointerB_C,
        &pointerE_C,
        &columns_C,
        &values_C);

    //
    valueResq = getValueFromExportSparse(r, c, pointerB_C, pointerE_C, columns_C, values_C);
    std::cout << "Find val row after set. Row " << r << " Colum " << c << " value " << valueResq << std::endl;
    //

    ///////////////////////////////////////////////////
    // Note: Should not be called before a set of the values (using mkl_sparse_s_set_value) in the sparse matrix
    status = mkl_sparse_set_mv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA, 1);
    mkl_sparse_optimize(csrA);
    ///////////////////////////////////////////////////

    status = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, csrA, descrA, vec, 0.0, &result2[0]);
    
    matrix[c + r * m] = valToSet;
    matrix[r + c * m] = valToSet;

    matmulVec(1.0, matrix, vec, 0.0, &result[0], m, m);

    for (int i = 0; i < m; ++i)
    {
        std::cout << " result after set with sgemv " << result[i] << " result after set with sparse_s_mv " << result2[i] << std::endl;
    }


}



