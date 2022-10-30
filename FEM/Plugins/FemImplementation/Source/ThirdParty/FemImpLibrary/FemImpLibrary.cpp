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

void matricesAdd(const float* A, const float* B, float* C, const int m, const int n)
{
    mkl_somatadd('r', 'n', 'n', m, n, 1.0f, A, n, 1.0f, B, n, C, n);
}

void matricesAdd(const double* A, const double* B, double* C, const int m, const int n)
{
    mkl_domatadd('r', 'n', 'n', m, n, 1.0, A, n, 1.0, B, n, C, n);
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
    float load(T t)
    {
        return 10000 * std::cos(T(1) * t);
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

    std::vector<T> storkm;
    std::vector<T> stormm;

    std::vector<T> u;
    std::vector<T> p;
    std::vector<T> xnew;


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

    std::vector<T> x0;
    std::vector<T> d1x0;
    std::vector<T> x1;
    std::vector<T> d2x0;
    std::vector<T> d1x1;
    std::vector<T> d2x1;

    T theta = T(0.5);
    T fm = T(0.005);
    T fk = T(0.012);

    T* loads2;

    int itt = 1;

    std::string   element;

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
    }

    void init(T* g_coord, int* g_num, int* in_nf, int in_nn)
    {
        nn = in_nn;
               
        //nprops = number of material properties
        const int nprops = 3;
        
        //np_types = number of diffent property types
        const int np_types = 1;
        
        nf.resize(nodof * nn, 0);
        
        std::copy(in_nf, in_nf + nodof * nn, std::begin(nf));
               
        neq = formnf(&nf[0], nodof, nn);
               
        //ndof = number of degree of freedom per element
        int ndof = nod * nodof;
        
        //prop = material property(e, v, gamma)
        std::vector<T> prop(nprops * np_types, 0.0f);
        
        prop[0] = T(1.0);
        prop[1] = T(0.3);
        prop[2] = T(1.0);
        
        //std::vector<int> g_g(ndof * nels, 0);
        g_g.resize(ndof * nels, 0);
        //kdiag.resize(neq, 0);
               
        //---------------------- loop the elements to find global arrays sizes---- -
        ///////////////
        T dtim = T(0.0);
        T c1 = (T(1.0) - theta) * dtim;
        T c2 = fk - c1;
        T c3 = fm + T(1.0) / (theta * dtim);
        T c4 = fk + theta * dtim;
        ///////////////

        std::vector<int> g(ndof, 0);
        
        for (int i = 0; i < nels; ++i)
        {
            int* num = &g_num[nod * i];
        
            num_to_g(num, &nf[0], &g[0], nod, nodof, nn);
        
            for (int j = 0; j < ndof; ++j)
            {
                g_g[i + j * nels] = g[j];
            }
        }
      
        sample(element, points, weights, ndim);

        gravlo.resize(neq, T(0.0));
        diag_precon.resize(neq + 1, T(0.0));

        //----element stiffnessand mass integration, storageand preconditioner-- -

        const T e = prop[0];
        const T v = prop[1];
        
        int* etype = new int[nels];
        std::fill(etype, etype + nels, 1);
        
        //km = element stiffness matrix
        std::vector<T> km(ndof * ndof, T(0.0));
        std::vector<T> jac  (ndim * ndim, T(0.0));
        std::vector<T> deriv(ndim * nod, T(0.0));
        std::vector<T> bee  (nst * ndof, T(0.0));
        //mm = element mass matrix;
        std::vector<T> mm(ndof * ndof, 0.0f);
        
        for (int i = 0; i < nels; ++i)
        {
            std::fill(std::begin(mm), std::end(mm), T(0.0));
        
            std::vector<T> dee(nst * nst, T(0.0));
            //std::fill(std::begin(dee), std::end(dee), T(0.0));
            deemat(&dee[0], nst, e, v);
        
            int* num = &g_num[nels * i];
        
            std::vector<T> coord(nod * ndim, T(0.0));
        
            for (int j = 0; j < nod; ++j)
            {
                for (int k = 0; k < ndim; ++k)
                {
                    coord[k + j * ndim] = g_coord[num[j] + k * nn];
                }
            }
        
            std::vector<int> g(ndof, 0);
            for (int j = 0; j < ndof; ++j)
            {
                g[j] = g_g[i + j * nels];
            }
        
            std::fill(std::begin(km), std::end(km), T(0.0));
            
            std::vector<T> mm(ndof * ndof, T(0.0));
            ////// for each point of integration, we have just one for 3d tets
            
            for (int j = 0; j < nip; ++j)
            {
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

                matricesAdd(&km[0], &kmTemp[0], &km[0], ndof, ndof);

                std::vector<T> ecm(ndof * ndof, T(0.0));
                std::vector<T> nt(ndof * nodof, T(0.0));
                std::vector<T> tn(nodof * ndof, T(0.0));

                ecmat(&ecm[0], &nt[0], &tn[0], &fun[0], ndof, nodof);

                for (int k = 0; k < ndof * ndof; ++k)
                {
                    mm[k] = mm[k] + ecm[k] * det * weights[j] * prop[2];
                }
            }

            //for (int j = nodof; j <= ndof; j += nodof)
            //{
            //    eld[j - 1] = fun[int(j / nodof) - 1] * det * weights;
            //}
        
            
            
            
            storkm.insert(storkm.end(), km.begin(), km.end());
            stormm.insert(stormm.end(), mm.begin(), mm.end());

            for (int k = 1; k < ndof; ++k)
            {
                diag_precon[g[k - 1]] = diag_precon[g[k - 1]] + mm[k + (k - 1) * ndof - 1] * c3 + km[k + (k - 1) * ndof - 1] * c4;
            }
        }

        for (int k = 1; k < diag_precon.size(); ++k)
        {
            diag_precon[k] = T(1.0) / diag_precon[k];
        }
        diag_precon[0] = T(0.0);
        //---------------------- - initial conditions and factorise equations--------

        // x0 = old displacements
        // x1 = new displacements
        // d1x0 = old velocity
        // d1x1 = new velocity
        // d2x0 = old acceleration
        // d2x1 = new acceleration
        // theta = time integration weighting parameter
        // fk = Rayleigh damping parameter on stiffness
        // fm = Rayleigh damping parameter on mass

        //
       

        x0.resize(neq + 1, T(0.0));
        d1x0.resize(neq + 1, T(0.0));
        //x1.resize(neq + 1, T(0.0));
        d2x0.resize(neq + 1, T(0.0));
        //d1x1.resize(neq + 1, T(0.0));
        //d2x1.resize(neq + 1, T(0.0));

               

        //T c1 = (T(1.0) - theta) * dtim;
        //T c2 = fk - c1;
        //T c3 = fm + T(1.0) / (theta * dtim);
        //T c4 = fk + theta * dtim;

        //addVectors(&mv[0], c3, &kv[0], c4, &f1[0], kdiag[neq - 1]);
        //sparin(&f1[0], &kdiag[0], neq);

        //!---------------------- - time stepping loop--------------------------------

        //T nstep = T(20);

        //std::cout << "time obj " << time << "      load  obj  " << load(time) << "     x   obj " << x0[nf[nres - 1]] << "     y  obj " << x0[nf[nres + nn - 1]] << "     z obj   " << x0[nf[nres + 2 * nn - 1]] << std::endl;
        
        //loads2 = new T[neq + 1];

        u.resize(neq + 1, T(0.0));
    }

    void update(T dtim, T* verticesBuffer)
    {
        ////number of loaded nodes
        //const int loaded_nodes = 2;
        //int node[loaded_nodes] = {};
        //
        ////val = applied nodal load weightings
        //std::vector<T> val(loaded_nodes * ndim, T(0.0));
        //
        //for (int i = 0; i < loaded_nodes; ++i)
        //{
        //    for (int j = 0; j < ndim; ++j)
        //    {
        //        val[i * loaded_nodes + j] = T(0.25);
        //    }
        //}
        //
        //// for debug purposes //
        ////nres = node number at witch time history is to be printed
        //int nres = 6;
        //node[0] = nres;
        //node[1] = 2;
        //int npri = 1;
        //////////////////////////
        //
        int ndof = nod * nodof;

        T c1 = (T(1.0) - theta) * dtim;
        T c2 = fk - c1;
        T c3 = fm + T(1.0) / (theta * dtim);
        T c4 = fk + theta * dtim;
        //
        //addVectors(&mv[0], c3, &kv[0], c4, &f1[0], kdiag[neq - 1]);
        //sparin(&f1[0], &kdiag[0], neq);
        //
        ////for (int i = 1; i <= nstep; ++i)
        ////{
        time = time + dtim;
        std::fill(loads2, loads2 + neq + 1, T(0.0));
        std::fill(u.begin(), u.end(), T(0.0));

        for (int i = 0; i < nels; ++i)
        {
            std::vector<T> x0Temp(ndof, 0);
            std::vector<T> d1x0Temp(ndof, 0);
            std::vector<T> uTemp(ndof, 0);

            for (int j = 0; j < ndof; ++j)
            {
                int g_index = g_g[i + j * nels];

                x0Temp[j] = x0[g_index];
                d1x0Temp[j] = d1x0[g_index];
            }

            T* km = &storkm[i * ndof * ndof];
            T* mm = &stormm[i * ndof * ndof];

            //u(g) = u(g) + MATMUL(km * c2 + mm * c3, x0(g)) + MATMUL(mm / theta, d1x0(g))

            //MATMUL(km * c2, X0(g))
            //MATMUL(mm * c3, x0(g))
            //MATMUL(mm / theta, d1x0(g))

            matmulVec(c2, km, &x0Temp[0], T(1.0), &uTemp[0], ndof, ndof);
            matmulVec(c3, mm, &x0Temp[0], T(1.0), &uTemp[0], ndof, ndof);
            matmulVec(T(1.0) / theta, mm, &d1x0Temp[0], T(1.0), &uTemp[0], ndof, ndof);
        }

        //
        //addVectors(&x0[0], c3, &d1x0[0], 1.0f / theta, &x1[0], neq + 1);
        //
        //float temporal = theta * dtim * load(time) + c1 * load(time - dtim);
        //
        //for (int j = 1; j <= loaded_nodes; ++j)
        //{
        //    for (int k = 0; k < ndim; ++k)
        //    {
        //        loads2[nf[nn * k + node[j - 1] - 1]] =
        //            val[(j - 1) * loaded_nodes + k] * temporal;
        //    }
        //}
        //
        //linmul_sky(&mv[0], &x1[0], &d1x1[0], &kdiag[0], neq);
        //
        ////d1x1=loads+d1x1
        //addVectors(T(1.0), loads2, &d1x1[0], neq + 1);
        //
        //copyVec(neq + 1, &x0[0], loads2);
        //scalVecProduct(neq + 1, c2, loads2);
        //
        //linmul_sky(&kv[0], loads2, &x1[0], &kdiag[0], neq);
        //
        ////x1=x1+d1x1
        //addVectors(T(1.0), &d1x1[0], &x1[0], neq + 1);
        //
        //spabac(&f1[0], &x1[0], &kdiag[0], neq);
        //
        //T a = T(1.0) / (theta * dtim);
        //T b = (T(1.0) - theta) / theta;
        //
        ////d1x1 = a * (x1 - x0) - b * d1x0;
        //addVectors(&x1[0], a, &x0[0], -a, &d1x1[0], neq + 1);
        //addVectors(-b, &d1x0[0], &d1x1[0], neq + 1);
        //
        ////d2x1 = a * (d1x1 - d1x0) - b * d2x0;
        //addVectors(&d1x1[0], a, &d1x0[0], -a, &d2x1[0], neq + 1);
        //addVectors(-b, &d2x0[0], &d2x1[0], neq + 1);
        //
        //copyVec(neq + 1, &x1[0], &x0[0]);
        //copyVec(neq + 1, &d1x1[0], &d1x0[0]);
        //copyVec(neq + 1, &d2x1[0], &d2x0[0]);
        //
        //if (itt / npri * npri == itt)
        //{
        ////  std::cout << "time obj2 " << time << "      load  obj  " << load(time) << "     x   obj " << x0[nf[nres - 1]] << "     y  obj " << x0[nf[nres + nn - 1]] << "     z obj   " << x0[nf[nres + 2 * nn - 1]] << std::endl;
        //}
        //++itt;
        //
        //if (verticesBuffer)
        //{
        //    for (int i = 0; i < nn; ++i)
        //    {
        //        verticesBuffer[i * 3] += (nf[i] > 0) ? x0[nf[i]] : T(0.0);
        //        verticesBuffer[i * 3 + 1] += (nf[i + nn] > 0) ? x0[nf[i + nn]] : T(0.0);
        //        verticesBuffer[i * 3 + 2] += (nf[i + 2 * nn] > 0) ? x0[nf[i + 2 * nn]] : T(0.0);
        //
        //        //if (i == 5)
        //        //{
        //        //    std::cout << x0[nf[i]] << " " << x0[nf[i + nn]] << " " << x0[nf[i + 2 * nn]] << std::endl;
        //        //}
        //    }
        //}
        ////}
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
void FEM_Factory<T>::update(int id, T dtim, T* verticesBuffer)
{
    if (id < femAlg.size())
    {
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



