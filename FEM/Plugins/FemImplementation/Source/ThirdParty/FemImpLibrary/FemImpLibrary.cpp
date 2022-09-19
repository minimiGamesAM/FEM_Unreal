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

    template<class T>
    void deemat(T* dee, int nbRowsDee, T e, T v)
    {
        T ih = nbRowsDee;
        T one(T(1.0)), two(T(2.0)), pt5 = T(0.5);
        T v1 = one - v;
        T c = e / ((one + v) * (one - two * v));

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

        for (int i = 0; i < nbRowsDee* nbRowsDee; ++i)
        {
            dee[i] *= e / (two * (one + v) * vv);
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

void addVectors(float* X, float a, float* Y, float b, float* R, float n)
{
    cblas_sscal(n, b, Y, 1);
    cblas_scopy(n, Y, 1, R, 1);
    cblas_saxpy(n, a, X, 1, R, 1);
}

void addVectors(double* X, double a, double* Y, double b, double* R, double n)
{
    cblas_dscal(n, b, Y, 1);
    cblas_dcopy(n, Y, 1, R, 1);
    cblas_daxpy(n, a, X, 1, R, 1);
}

void addVectors(float a, float* X, float* Y, float n)
{
    cblas_saxpy(n, a, X, 1, Y, 1);
}

void addVectors(double a, double* X, double* Y, double n)
{
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

// C = alpha A * B + beta * C
void matmul(float* A, float* B, float* C, int m, int k, int n)
{
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);
}

void matmul(double* A, double* B, double* C, int m, int k, int n)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, k, B, n, 0, C, n);
}

void matmulTransA(float* A, float* B, float* C, int m, int k, int n)
{
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, B, n, 0, C, n);
}

void matmulTransA(double* A, double* B, double* C, int m, int k, int n)
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
        return std::cos(T(0.3) * t);
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
    const int nod = 4;

    //nf = nodal freedom array(nodof rows and nn colums)
    int* nf;
    
    T fun[4] = { T(0.25), T(0.25), T(0.25), T(0.25) };

    T der[3 * 4] = 
    { //[ndim * nod] = {
            T(1.0), T(0.0), T(0.0), T(-1.0),
            T(0.0), T(1.0), T(0.0), T(-1.0),
            T(0.0), T(0.0), T(1.0), T(-1.0)
    };

    // skyline profile
    std::vector<int> kdiag;

    //kv = global stiffness matrix
    std::vector<T> kv;
   
    //global consisten mass
    std::vector<T> mv;

    //left hand side matrix (stored as a skyline) 
    std::vector<T> f1;

    //gravlo = global gravity loading vector 
    std::vector<T> gravlo;

public:
    Fem_Algoritm(int ndim, int nodof, int nels)
        : ndim(ndim), nodof(nodof), nn(0), nels(nels), neq(0)
    {
    }

    void init(T* g_coord, int* g_num, int* in_nf, int in_nn)
    {
        nn = in_nn;
        
        //nip = number of intregation points per element
        const int nip = 1;
        
        //nprops = number of material properties
        const int nprops = 3;
        
        //np_types = number of diffent property types
        const int np_types = 1;
        
        nf = new int[nodof * nn];
        
        std::copy(in_nf, in_nf + nodof * nn, nf);
               
        neq = formnf(nf, nodof, nn);
               
        //ndof = number of degree of freedom per element
        int ndof = nod * nodof;
        
        //prop = material property(e, v, gamma)
        std::vector<T> prop(nprops * np_types, 0.0f);
        
        prop[0] = T(100.0);
        prop[1] = T(0.3);
        prop[2] = T(1.0);
        
        std::vector<int> g_g(ndof * nels, 0);
               
        kdiag.resize(neq, 0);
               
        //---------------------- loop the elements to find global arrays sizes---- -
        
        std::vector<int> g(ndof, 0);
        
        for (int i = 0; i < nels; ++i)
        {
            int* num = &g_num[4 * i];
        
            num_to_g(num, nf, &g[0], nod, nodof, nn);
        
            for (int j = 0; j < ndof; ++j)
            {
                g_g[i + j * nels] = g[j];
            }
        
            fkdiag(&kdiag[0], &g[0], ndof);
        }
        
        for (int i = 1; i < neq; ++i)
        {
            kdiag[i] = kdiag[i] + kdiag[i - 1];
        }
        
        kv.resize(kdiag[neq - 1], T(0.0));
        mv.resize(kdiag[neq - 1], T(0.0));
        f1.resize(kdiag[neq - 1], T(0.0));
               
        gravlo.resize(neq, T(0.0));
         
        //----------------------- element stiffness integration and assembly--------
        
        //call sample, but for tet there is just one point
        T points[] = { T(0.25), T(0.25), T(0.25) };
        T weights = T(1.0) / T(6.0);
        
        //nst = number of stress / strain terms
        const int nst = 6;
        const T e = T(100.0);
        const T v = T(0.3);
        
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
        
            T dee[nst * nst] = {};
            std::fill(std::begin(dee), std::end(dee), T(0.0));
            deemat(dee, nst, e, v);
        
            int* num = &g_num[4 * i];
        
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
            // calculate jac
            matmul(der, &coord[0], &jac[0], ndim, nod, ndim);
            T det = invert(&jac[0], ndim);
        
            // calculate the derivative in x, y, z
            matmul(&jac[0], der, &deriv[0], ndim, ndim, nod);
        
            std::fill(std::begin(bee), std::end(bee), T(0.0));
            beemat(&bee[0], nst, ndof, &deriv[0], nod);
        
            std::vector<T> temporal(nst * ndof, T(0.0));
                    
            matmulTransA(&bee[0], dee, &temporal[0], ndof, nst, nst);
            matmul(&temporal[0], &bee[0], &km[0], ndof, nst, ndof);
        
            std::for_each(std::begin(km), std::end(km), [&](T& v) { v *= det * weights; });
        
            //for (int j = nodof; j <= ndof; j += nodof)
            //{
            //    eld[j - 1] = fun[int(j / nodof) - 1] * det * weights;
            //}
        
            std::vector<T> ecm(ndof* ndof, T(0.0));
            std::vector<T> nt (ndof * nodof, T(0.0));
            std::vector<T> tn (nodof * ndof, T(0.0));
                
            ecmat(&ecm[0], &nt[0], &tn[0], fun, ndof, nodof);
        
            for (int j = 0; j < ndof * ndof; ++j)
            {
                mm[j] = mm[j] + ecm[j] * det * weights * prop[2];
            }
        
            ////// end for each integration point
        
            fsparv(&kv[0], &km[0], &g[0], &kdiag[0], ndof);
            fsparv(&mv[0], &mm[0], &g[0], &kdiag[0], ndof);
        }
    }

    void update()
    {
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
        T dtim  = T(0.2);
        T theta = T(0.5);
        T fm    = T(0.005);
        T fk    = T(0.272);

        std::vector<T> x0   (neq + 1, T(0.0));
        std::vector<T> d1x0 (neq + 1, T(0.0));
        std::vector<T> x1   (neq + 1, T(0.0));
        std::vector<T> d2x0 (neq + 1, T(0.0));
        std::vector<T> d1x1 (neq + 1, T(0.0));
        std::vector<T> d2x1 (neq + 1, T(0.0));

        //number of loaded nodes
        const int loaded_nodes = 2;
        int node[loaded_nodes] = {};

        //val = applied nodal load weightings
        std::vector<T> val(loaded_nodes * ndim, T(0.0));

        for (int i = 0; i < loaded_nodes; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                val[i * loaded_nodes + j] = T(0.25);
            }
        }

        // for debug purposes //
        //nres = node number at witch time history is to be printed
        int nres = 6;
        node[0] = nres;
        node[1] = 2;
        ////////////////////////

        T c1 = (T(1.0) - theta) * dtim;
        T c2 = fk - c1;
        T c3 = fm + T(1.0) / (theta * dtim);
        T c4 = fk + theta * dtim;

        addVectors(&mv[0], c3, &kv[0], c4, &f1[0], kdiag[neq - 1]);
        sparin(&f1[0], &kdiag[0], neq);

        //!---------------------- - time stepping loop--------------------------------

        T time = T(0);
        T nstep = T(20);

        std::cout << "time obj " << time << "      load  obj  " << load(time) << "     x   obj " << x0[nf[nres - 1]] << "     y  obj " << x0[nf[nres + nn - 1]] << "     z obj   " << x0[nf[nres + 2 * nn - 1]] << std::endl;
       
        int npri = 1;

        T* loads2 = new T[neq + 1];

        for (int i = 1; i <= nstep; ++i)
        {
            time = time + dtim;
            std::fill(loads2, loads2 + neq + 1, T(0.0));

            addVectors(&x0[0], c3, &d1x0[0], 1.0f / theta, &x1[0], neq + 1);

            float temporal = theta * dtim * load(time) + c1 * load(time - dtim);

            for (int j = 1; j <= loaded_nodes; ++j)
            {
                for (int k = 0; k < ndim; ++k)
                {
                    loads2[nf[nn * k + node[j - 1] - 1]] =
                        val[(j - 1) * loaded_nodes + k] * temporal;
                }
            }

            linmul_sky(&mv[0], &x1[0], &d1x1[0], &kdiag[0], neq);

            //d1x1=loads+d1x1
            addVectors(T(1.0), loads2, &d1x1[0], neq + 1);

            copyVec(neq + 1, &x0[0], loads2);
            scalVecProduct(neq + 1, c2, loads2);

            linmul_sky(&kv[0], loads2, &x1[0], &kdiag[0], neq);

            //x1=x1+d1x1
            addVectors(T(1.0), &d1x1[0], &x1[0], neq + 1);

            spabac(&f1[0], &x1[0], &kdiag[0], neq);

            T a = T(1.0) / (theta * dtim);
            T b = (T(1.0) - theta) / theta;

            //d1x1 = a * (x1 - x0) - b * d1x0;
            addVectors(&x1[0], a, &x0[0], -a, &d1x1[0], neq + 1);
            addVectors(-b, &d1x0[0], &d1x1[0], neq + 1);

            //d2x1 = a * (d1x1 - d1x0) - b * d2x0;
            addVectors(&d1x1[0], a, &d1x0[0], -a, &d2x1[0], neq + 1);
            addVectors(-b, &d2x0[0], &d2x1[0], neq + 1);

            copyVec(neq + 1, &x1[0], &x0[0]);
            copyVec(neq + 1, &d1x1[0], &d1x0[0]);
            copyVec(neq + 1, &d2x1[0], &d2x0[0]);

            if (i / npri * npri == i)
            {
                std::cout << "time obj2 " << time << "      load  obj  " << load(time) << "     x   obj " << x0[nf[nres - 1]] << "     y  obj " << x0[nf[nres + nn - 1]] << "     z obj   " << x0[nf[nres + 2 * nn - 1]] << std::endl;
            }
        }
    }
};

template<class T>
void FEM_Factory<T>::create(int ndim, int nodof, int nels)
{
    femAlg = new Fem_Algoritm<T>(ndim, nodof, nels);
}

template<class T>
void FEM_Factory<T>::init(T* g_coord, int* g_num, int* in_nf, int in_nn)
{
    femAlg->init(g_coord, g_num, in_nf, in_nn);
}

template<class T>
void FEM_Factory<T>::update()
{
    femAlg->update();
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



