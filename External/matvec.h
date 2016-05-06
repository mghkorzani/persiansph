/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_MATVEC_H
#define MECHSYS_MATVEC_H

// Std Lib
#include <iostream>
#include <sstream>   // for istringstream, ostringstream
#include <cmath>     // for sqrt, pow
#include <algorithm> // for min, max

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// Tensors
#ifdef HAS_TENSORS
  #include <tensors/operators.h>
  #include <tensors/tensor1.h>
  #include <tensors/tensor2.h>
  #include <tensors/tensor3.h>
  #include <tensors/tensor4.h>
  typedef TensorsLib::Tensor1<double,3> Ten1_t;
  typedef TensorsLib::Tensor2<double,3> Ten2_t;
  typedef TensorsLib::Tensor3<double,3> Ten3_t;
  typedef TensorsLib::Tensor4<double,3> Ten4_t;
#endif

// MechSys
#include "fatal.h"
#include "util.h"
#include "array.h"

// LAPACK
extern "C"
{
    // DGETRF - compute an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges
    void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);

    // DGETRI - compute the inverse of a matrix using the LU factorization computed by DGETRF
    void dgetri_(const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);

    // DGESVD - computes the singular value decomposition of a real M-by-N matrix A, optionally computing the left and/or right singular vectors
    void dgesvd_(const char* jobu, const char* jobvt, const int* M, const int* N, double* A, const int* lda, double* S, double* U, const int* ldu, double* VT, const int* ldvt, double* work,const int* lwork, const int* info);

    // DGESV - double-general-solver
    void dgesv_(int *Np, int *NRHSp, double *A, int *LDAp, int *IPIVp, double *B, int *LDBp, int *INFOp);
}


/////////////////////////////////////////////////////////////////////////////////////////// General MatVec /////


#ifdef USE_MTL4

// Boost/MTL4
#include <boost/numeric/mtl/mtl.hpp>

/** Dense matrix (general). */
typedef mtl::dense2D<double, mtl::matrix::parameters<mtl::tag::col_major> > Mat_t;

/** Dense vector (general). */
typedef mtl::dense_vector<double> Vec_t;

#else

#include "../External/vector.h"
#include "../External/matrix.h"
#include "../External/laexpr.h"

typedef LinAlg::Vector<double> Vec_t;
typedef LinAlg::Matrix<double> Mat_t;

inline double dot         (Vec_t const & V, Vec_t const & W) { return LinAlg::Dot(V,W); }
inline size_t size        (Vec_t const & V) { return V.Size(); }
inline void   set_to_zero (Vec_t       & V) { V.SetValues(0.0); }
inline void   set_to_zero (Mat_t       & M) { M.SetValues(0.0); }
inline size_t num_rows    (Vec_t const & V) { return V.Size(); }
inline size_t num_rows    (Mat_t const & M) { return M.Rows(); }
inline size_t num_cols    (Mat_t const & M) { return M.Cols(); }

#endif

/** Print vector. */
inline String PrintVector (Vec_t const & V, char const * Fmt="%13g", Array<long> const * SkipR=NULL, double Tol=1.0e-13)
{
    int m = size(V);
    String lin;
    for (int i=0; i<m; ++i)
    {
        bool skip_row = false;
        if (SkipR!=NULL) skip_row = (SkipR->Find(i)<0 ? false : true);
        if (!skip_row)
        {
            double val = (fabs(V(i))<Tol ? 0.0 : V(i));
            String buf;  buf.Printf(Fmt,val);
            lin.append(buf);
        }
    }
    lin.append("\n");
    return lin;
}

/** Print matrix. */
inline String PrintMatrix (Mat_t const & M, char const * Fmt="%13g", Array<long> const * SkipRC=NULL, double Tol=1.0e-13, bool NumPy=false)
{
    int m = M.num_rows();
    int n = M.num_cols();
    String lin;
    if (NumPy) lin.append("matrix([[");
    for (int i=0; i<m; ++i)
    {
        bool skip_row = false;
        if (SkipRC!=NULL) skip_row = (SkipRC->Find(i)<0 ? false : true);
        if (!skip_row)
        {
            if (NumPy && i!=0) lin.append("        [");
            for (int j=0; j<n; ++j)
            {
                bool skip_col = false;
                if (SkipRC!=NULL) skip_col = (SkipRC->Find(j)<0 ? false : true);
                if (!skip_col)
                {
                    double val = (fabs(M(i,j))<Tol ? 0.0 : M(i,j));
                    String buf;  buf.Printf(Fmt,val);
                    lin.append(buf);
                    if (NumPy && j!=n-1) lin.append(",");
                }
            }
            if (NumPy && i!=m-1) lin.append("],");
            if (NumPy && i==m-1) lin.append("]])");
            else lin.append("\n");
        }
    }
    return lin;
}

/** Write SMatrix for vismatrix. */
inline void WriteSMAT (Mat_t const & M, char const * FileKey, double Tol=1.0e-14)
{
    // find the number of really non-zero values
    size_t m  = num_rows(M);
    size_t n  = num_cols(M);
    size_t nz = 0;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
    {
        if (fabs(M(i,j))>Tol) nz++;
    }

    // output
    std::ostringstream oss;
    char buf[256];
    sprintf(buf, "%zd  %zd  %zd\n", m, n, nz);
    oss << buf;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
    {
        if (fabs(M(i,j))>Tol)
        {
            sprintf(buf, "  %zd  %zd  %23.15e\n", i, j, M(i,j));
            oss << buf;
        }
    }

    // write to file
    String fn(FileKey);  fn.append(".smat");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
    printf("File <%s%s%s> written\n",TERM_CLR_BLUE_H, fn.CStr(), TERM_RST);
}

/** Compare two vectors. */
inline double CompareVectors (Vec_t const & A, Vec_t const & B)
{
    size_t m = size(A);
    if (m!=size(B)) throw new Fatal("matvec.h:CompareVectors: vectors A_%d and B_%d must have the same size",m,size(B));
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
        error += fabs(A(i)-B(i));
    return error;
}

/** Compare two matrices. */
inline double CompareMatrices (Mat_t const & A, Mat_t const & B)
{
    size_t m = A.num_rows();
    size_t n = A.num_cols();
    if ((m!=B.num_rows()) || (n!=B.num_cols())) throw new Fatal("matvec.h:CompareMatrices: matrices A_%dx%d and B_%dx%d must have the same number of rows and columns",m,n,B.num_rows(),B.num_cols());
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
        error += fabs(A(i,j)-B(i,j));
    return error;
}

/** Check if matrix is diagonal. */
inline double CheckDiagonal (Mat_t const & M, bool CheckUnitDiag=false)
{
    size_t m = M.num_rows();
    size_t n = M.num_cols();
    double error = 0.0;
    for (size_t i=0; i<m; ++i)
    for (size_t j=0; j<n; ++j)
    {
        if ((i==j) && CheckUnitDiag) error += fabs(M(i,j)-1.0);
        if  (i!=j)                   error += fabs(M(i,j));
    }
    return error;
}

/** Determinant. */
inline double Det (Mat_t const & M)
{
    int m = M.num_rows();
    int n = M.num_cols();
    if (m==1)
    {
        double res = 0;
        for (int i=0; i<n; i++) res += M(0,i)*M(0,i);
        return sqrt(res);
    }
    else if (m==2 && n==2)
    {
        return M(0,0)*M(1,1) - M(1,0)*M(0,1);
    }
    /*
    else if (m==2 && n==3)
    {
        double d1 = M(0,0)*M(1,1) - M(0,1)*M(1,0);
        double d2 = M(0,1)*M(1,2) - M(0,2)*M(1,1);
        double d3 = M(0,2)*M(1,0) - M(0,0)*M(1,2);
        return sqrt(d1*d1 + d2*d2 + d3*d3);
    }
    */
    else if (m==3 && n==3)
    {
        return  M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
              - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
              + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    }
    else if (m==n)
    {
        // factorization
        int   info = 0;
        int * ipiv = new int [m];
        Mat_t Mcpy(M);
        dgetrf_(&m,         // M
                &m,         // N
                Mcpy.data,  // double * A
                &m,         // LDA
                ipiv,       // Pivot indices
                &info);     // INFO
        if (info!=0) throw new Fatal ("matvec.h::Det: LAPACK: LU factorization failed");

        // determinant
        double det = 1.0;
        for (int i=0; i<m; ++i)
        {
            if (ipiv[i]!=(i+1)) det = -det * Mcpy(i,i);
            else                det =  det * Mcpy(i,i);
        }

        // end
        delete [] ipiv;
        return det;
    }
    else throw new Fatal("matvec.h:Det: Method is not implemented for (%d x %d) matrices yet",m,n);
}

/** Identity. */
inline void Identity (size_t NRows, Mat_t & I)
{
    I.change_dim (NRows,NRows);
    for (size_t i=0; i<NRows; ++i)
    for (size_t j=0; j<NRows; ++j)
    {
        if (i==j) I(i,j) = 1.0;
        else      I(i,j) = 0.0;
    }
}

/** Singular value decomposition. M = U_mxm * D_mxn * Vt_nxn   */
inline void Svd (Mat_t const & M, Mat_t & U, Vec_t & S, Mat_t & Vt)
{
    int  info   = 0;
    char job    = 'A';
    int  m      = M.num_rows();
    int  n      = M.num_cols();
    int  min_mn = (m<n ? m : n);
    int  max_mn = (m>n ? m : n);
    int  lwork  = 2.0*std::max(3*min_mn+max_mn, 5*min_mn);

    U. change_dim (m, m);
    Vt.change_dim (n, n); // trans(V)
    S. change_dim (min_mn);

    double * work  = new double [lwork]; // Work

    // decomposition
    Mat_t tmp(M);
    dgesvd_(&job,      // JOBU
            &job,      // JOBVT
            &m,        // M
            &n,        // N
            tmp.data,  // A
            &m,        // LDA
            S.data,    // S
            U.data,    // U
            &m,        // LDU
            Vt.data,   // VT
            &n,        // LDVT
            work,      // WORK
            &lwork,    // LWORK
            &info);    // INFO
    if (info!=0) throw new Fatal ("matvec::Svd: LAPACK: Decomposition failed");

    delete [] work;
}

/** Inverse. */
inline void Inv (Mat_t const & M, Mat_t & Mi, double Tol=1.0e-10)
{
    int m = M.num_rows();
    int n = M.num_cols();
    Mi.change_dim(m,n);
    if (m==2 && n==2)
    {
        double det = Det(M);
        if (fabs(det)<Tol) throw new Fatal("matvec.h:Inv: Cannot calculate inverse due to zero determinant. det = %g",det);

        Mi(0,0) =  M(1,1) / det;
        Mi(0,1) = -M(0,1) / det;

        Mi(1,0) = -M(1,0) / det;
        Mi(1,1) =  M(0,0) / det;
    }
    else if (m==3 && n==3)
    {
        double det = Det(M);
        if (fabs(det)<Tol) throw new Fatal("matvec.h:Inv: Cannot calculate inverse due to zero determinant. det = %g",det);

        Mi(0,0) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det;
        Mi(0,1) = (M(0,2)*M(2,1) - M(0,1)*M(2,2)) / det;
        Mi(0,2) = (M(0,1)*M(1,2) - M(0,2)*M(1,1)) / det;

        Mi(1,0) = (M(1,2)*M(2,0) - M(1,0)*M(2,2)) / det;
        Mi(1,1) = (M(0,0)*M(2,2) - M(0,2)*M(2,0)) / det;
        Mi(1,2) = (M(0,2)*M(1,0) - M(0,0)*M(1,2)) / det;

        Mi(2,0) = (M(1,0)*M(2,1) - M(1,1)*M(2,0)) / det;
        Mi(2,1) = (M(0,1)*M(2,0) - M(0,0)*M(2,1)) / det;
        Mi(2,2) = (M(0,0)*M(1,1) - M(0,1)*M(1,0)) / det;
    }
    else if (m==n) // square
    {
        int   info = 0;
        int * ipiv = new int [m];

        // factorization
        Mi = M;
        dgetrf_(&m,       // M
                &m,       // N
                Mi.data,  // double * A
                &m,       // LDA
                ipiv,     // Pivot indices
                &info);   // INFO
        if (info!=0) throw new Fatal ("matvec.h::Inv: LAPACK: LU factorization failed");

        int      NB    = 4;                  // Optimal blocksize ?
        int      lwork = m*NB;               // Dimension of work >= max(1,m), optimal=m*NB
        double * work  = new double [lwork]; // Work

        // inversion
        dgetri_(&m,       // N
                Mi.data,  // double * A
                &m,       // LDA
                ipiv,     // Pivot indices
                work,     // work
                &lwork,   // dimension of work
                &info);   // INFO
        if (info!=0) throw new Fatal ("matvec::Inv: LAPACK: Inversion failed");

        delete [] ipiv;
        delete [] work;
    }
    else // generalized (pseudo) inverse
    {
        Mat_t U; Vec_t S; Mat_t Vt;
        Svd(M, U, S, Vt);
        Mat_t Di(m,n);
        set_to_zero(Di);
        for (size_t i=0; i<size(S); ++i)
        {
            if (S(i)>Tol) Di(i,i) = 1.0/S(i);
        }
        Mat_t tmp(U*Di*Vt);
        Mi.change_dim(n,m);
        Mi = trans(tmp);
    }
}

/** Linear Solver. {X} = [M]^{-1}{X} (M is lost) (X initially has the contents of the right-hand side) */
inline void Sol (Mat_t & M, Vec_t & X)
{
    int  m = M.num_rows();
    int  n = M.num_cols();
    int mv = size(X);
    if (m!=n)  throw new Fatal("Sol: Matrix must be square");
    if (m!=mv) throw new Fatal("Sol: Vector X must have the same number of rows of matrix M");

    int   info = 0;
    int   nrhs = 1; // vector X has 1 column
    int * ipiv = new int [m];
    dgesv_(&m,      // A(m,m)
           &nrhs,   // {X}(m,1) (RHS: Right Hand Side) 
           M.data,  // double * A
           &m,      // LDA
           ipiv,    // Pivot Indices
           X.data,  // double * Y
           &m,      // LDY
           &info);  // info
    delete [] ipiv;

    if (info!=0) 
    {
        throw new Fatal ("Sol: Linear solver (DGESV) failed (singular matrix?)");
    }
}

/** Linear Solver. {X} = [M]^{-1}{B}  */
inline void Sol (Mat_t const & M, Vec_t const & B, Vec_t & X)
{
    Mat_t m(M);
    X = B;
    Sol (m, X);
}

/** Norm. */
inline double Norm (Vec_t const & V)
{
    return sqrt(dot(V,V));
}

/** Dyadic product. */
inline void Dyad (Vec_t const & A, Vec_t const & B, Mat_t & M)
{
    M.change_dim(size(A),size(B));
    for (size_t i=0; i<size(A); ++i)
    for (size_t j=0; j<size(B); ++j)
        M(i,j) = A(i) * B(j);
}

/** Left multiplication. {B} = {A}*[M]. NOTE: this is not efficient for large matrices.  */
inline void Mult (Vec_t const & A, Mat_t const & M, Vec_t & B)
{
    B.change_dim (M.num_cols());
    set_to_zero  (B);
    for (size_t i=0; i<M.num_rows(); ++i)
    for (size_t j=0; j<M.num_cols(); ++j)
        B(j) += A(i)*M(i,j);
}

/** Create column matrix from vector. */
inline void Vec2ColMat (Vec_t const & V, Mat_t & M)
{
    M.change_dim (size(V),1);
    for (size_t i=0; i<size(V); ++i) M(i,0) = V(i);
}

/** Eigenvalues of symmetric matrix. NOTE: This function changes the matrix M. */
inline void Eig (Mat_t & M, Vec_t & L, Mat_t * Q=NULL, bool Qtrans=false)
{
    // calculate
    size_t nrow = num_rows (M);
    size_t ncol = num_cols (M);
    if (nrow!=ncol) throw new Fatal("Eig: Matrix (%zdx%zd) must be square", nrow,ncol);

    // eigenvalues
    L.change_dim (nrow);
    gsl_vector eval = {nrow, 1, L.data, NULL, 0}; // size, stride, data, block, owner

    // solve
    gsl_matrix_view m = gsl_matrix_view_array (M.data, nrow, nrow);
    if (Q==NULL)
    {
        gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (nrow);
        gsl_eigen_symm (&m.matrix, &eval, w);
        gsl_eigen_symm_free (w);
    }
    else
    {
        Q->change_dim (nrow, nrow);
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (nrow);
        if (Qtrans)
        {
            gsl_matrix_view q = gsl_matrix_view_array (Q->data, nrow, nrow);
            gsl_eigen_symmv (&m.matrix, &eval, &q.matrix, w);
        }
        else
        {
            gsl_matrix * q = gsl_matrix_alloc (nrow, nrow);
            gsl_eigen_symmv (&m.matrix, &eval, q, w);
            for (size_t j=0; j<ncol; ++j)
            {
                gsl_vector_view evec_j = gsl_matrix_column (q, j);
                for (size_t i=0; i<nrow; ++i) (*Q)(i,j) = gsl_vector_get (&evec_j.vector, i);
            }
            gsl_matrix_free (q);
        }
        gsl_eigen_symmv_free (w);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////// Tiny MatVec /////


/** 3x3 Matrix. */
typedef blitz::TinyMatrix<double,3,3> Mat3_t;

/** 3x1 Vector. */
typedef blitz::TinyVector<double,3> Vec3_t;
typedef blitz::TinyVector<size_t,3> iVec3_t;
typedef blitz::TinyVector<bool,3>   bVec3_t;

Mat3_t operator * ( double a, const Mat3_t & A)
{
    Mat3_t M;
    M(0,0)=A(0,0)*a;  M(0,1)=A(0,1)*a;  M(0,2)=A(0,2)*a;
    M(1,0)=A(1,0)*a;  M(1,1)=A(1,1)*a;  M(1,2)=A(1,2)*a;
    M(2,0)=A(2,0)*a;  M(2,1)=A(2,1)*a;  M(2,2)=A(2,2)*a;
    return M;
}

Mat3_t operator + (const Mat3_t & A, const Mat3_t & B)
{
    Mat3_t M;
    M(0,0)=A(0,0)+B(0,0);  M(0,1)=A(0,1)+B(0,1);  M(0,2)=A(0,2)+B(0,2);
    M(1,0)=A(1,0)+B(1,0);  M(1,1)=A(1,1)+B(1,1);  M(1,2)=A(1,2)+B(1,2);
    M(2,0)=A(2,0)+B(2,0);  M(2,1)=A(2,1)+B(2,1);  M(2,2)=A(2,2)+B(2,2);
    return M;
}

Mat3_t operator - (const Mat3_t & A, const Mat3_t & B)
{
    Mat3_t M;
    M(0,0)=A(0,0)-B(0,0);  M(0,1)=A(0,1)-B(0,1);  M(0,2)=A(0,2)-B(0,2);
    M(1,0)=A(1,0)-B(1,0);  M(1,1)=A(1,1)-B(1,1);  M(1,2)=A(1,2)-B(1,2);
    M(2,0)=A(2,0)-B(2,0);  M(2,1)=A(2,1)-B(2,1);  M(2,2)=A(2,2)-B(2,2);
    return M;
}

/** Print vector. */
inline String PrintVector (Vec3_t const & V, char const * Fmt="%13g", double Tol=1.0e-13)
{
    int m = 3;
    String lin;
    for (int i=0; i<m; ++i)
    {
        double val = (fabs(V(i))<Tol ? 0.0 : V(i));
        String buf;  buf.Printf(Fmt,val);
        lin.append(buf);
    }
    lin.append("\n");
    return lin;
}

/** Print matrix. */
inline String PrintMatrix (Mat3_t const & M, char const * Fmt="%13g", double Tol=1.0e-13)
{
    int m = 3;
    int n = 3;
    String lin;
    for (int i=0; i<m; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            double val = (fabs(M(i,j))<Tol ? 0.0 : M(i,j));
            String buf;  buf.Printf(Fmt,val);
            lin.append(buf);
        }
        lin.append("\n");
    }
    return lin;
}

/** Compare two vectors. */
inline double CompareVectors (Vec3_t const & A, Vec3_t const & B)
{
    double error = 0.0;
    for (size_t i=0; i<3; ++i)
        error += fabs(A(i)-B(i));
    return error;
}

/** Compare two matrices. */
inline double CompareMatrices (Mat3_t const & A, Mat3_t const & B)
{
    double error = 0.0;
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
        error += fabs(A(i,j)-B(i,j));
    return error;
}

/** Transpose.*/
inline void Trans (Mat3_t const & M, Mat3_t & Mt)
{
    Mt(0,0)=M(0,0);   Mt(0,1)=M(1,0);   Mt(0,2)=M(2,0);
    Mt(1,0)=M(0,1);   Mt(1,1)=M(1,1);   Mt(1,2)=M(2,1);
    Mt(2,0)=M(0,2);   Mt(2,1)=M(1,2);   Mt(2,2)=M(2,2);
}

/** Determinant.*/
inline double Det (Mat3_t const & M)
{
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    return det;
}

/** Identity. */
inline void Identity (Mat3_t & I)
{
    I = 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
}

/** Inverse.*/
inline void Inv (Mat3_t const & M, Mat3_t & Mi, double Tol=1.0e-10)
{
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));

    if (fabs(det)<Tol)
    {
        std::ostringstream oss;
        oss << PrintMatrix(M);
        throw new Fatal("matvec.h::Inv: 3x3 matrix inversion failed with null (%g) determinat. M =\n%s",Tol,oss.str().c_str());
    }

    Mi(0,0)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))/det;  Mi(0,1)=(M(0,2)*M(2,1)-M(0,1)*M(2,2))/det;  Mi(0,2)=(M(0,1)*M(1,2)-M(0,2)*M(1,1))/det;
    Mi(1,0)=(M(1,2)*M(2,0)-M(1,0)*M(2,2))/det;  Mi(1,1)=(M(0,0)*M(2,2)-M(0,2)*M(2,0))/det;  Mi(1,2)=(M(0,2)*M(1,0)-M(0,0)*M(1,2))/det;
    Mi(2,0)=(M(1,0)*M(2,1)-M(1,1)*M(2,0))/det;  Mi(2,1)=(M(0,1)*M(2,0)-M(0,0)*M(2,1))/det;  Mi(2,2)=(M(0,0)*M(1,1)-M(0,1)*M(1,0))/det;
}

/** Linear Solver. {X} = [M]^{-1}{B}  */
inline void Sol (Mat3_t const & M, Vec3_t const & B, Vec3_t & X, double Tol=1.0e-10)
{
    // determinant
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    if (fabs(det)<Tol) throw new Fatal("matvec.h:Sol: Cannot calculate inverse due to zero determinant. det = %g",det);

    // inverse matrix
    Mat3_t Mi;
    Mi(0,0) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det;
    Mi(0,1) = (M(0,2)*M(2,1) - M(0,1)*M(2,2)) / det;
    Mi(0,2) = (M(0,1)*M(1,2) - M(0,2)*M(1,1)) / det;

    Mi(1,0) = (M(1,2)*M(2,0) - M(1,0)*M(2,2)) / det;
    Mi(1,1) = (M(0,0)*M(2,2) - M(0,2)*M(2,0)) / det;
    Mi(1,2) = (M(0,2)*M(1,0) - M(0,0)*M(1,2)) / det;

    Mi(2,0) = (M(1,0)*M(2,1) - M(1,1)*M(2,0)) / det;
    Mi(2,1) = (M(0,1)*M(2,0) - M(0,0)*M(2,1)) / det;
    Mi(2,2) = (M(0,0)*M(1,1) - M(0,1)*M(1,0)) / det;

    // solve system
    X(0) = Mi(0,0)*B(0) + Mi(0,1)*B(1) + Mi(0,2)*B(2);
    X(1) = Mi(1,0)*B(0) + Mi(1,1)*B(1) + Mi(1,2)*B(2);
    X(2) = Mi(2,0)*B(0) + Mi(2,1)*B(1) + Mi(2,2)*B(2);
}

/** Alternative Solver that does not throw a fatal error when the matrix is singular */

inline bool SolAlt (Mat3_t const & M, Vec3_t const & B, Vec3_t & X, double Tol=1.0e-10)
{
    // determinant
    double det =   M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
                 - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
                 + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
    if (fabs(det)<Tol) return false;
    
    // inverse matrix
    Mat3_t Mi;
    Mi(0,0) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det;
    Mi(0,1) = (M(0,2)*M(2,1) - M(0,1)*M(2,2)) / det;
    Mi(0,2) = (M(0,1)*M(1,2) - M(0,2)*M(1,1)) / det;

    Mi(1,0) = (M(1,2)*M(2,0) - M(1,0)*M(2,2)) / det;
    Mi(1,1) = (M(0,0)*M(2,2) - M(0,2)*M(2,0)) / det;
    Mi(1,2) = (M(0,2)*M(1,0) - M(0,0)*M(1,2)) / det;

    Mi(2,0) = (M(1,0)*M(2,1) - M(1,1)*M(2,0)) / det;
    Mi(2,1) = (M(0,1)*M(2,0) - M(0,0)*M(2,1)) / det;
    Mi(2,2) = (M(0,0)*M(1,1) - M(0,1)*M(1,0)) / det;

    // solve system
    X(0) = Mi(0,0)*B(0) + Mi(0,1)*B(1) + Mi(0,2)*B(2);
    X(1) = Mi(1,0)*B(0) + Mi(1,1)*B(1) + Mi(1,2)*B(2);
    X(2) = Mi(2,0)*B(0) + Mi(2,1)*B(1) + Mi(2,2)*B(2);

    return true;
}


/** Eigenvalues. NOTE: This function changes the matrix M. */
inline void Eig (Mat3_t & M, Vec3_t & L)
{
    // calculate
    gsl_matrix_view m = gsl_matrix_view_array (M.data(), 3, 3);
    gsl_vector * eval = gsl_vector_alloc      (3);
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (3);
    gsl_eigen_symm (&m.matrix, eval, w);

    // eigenvalues
    L = gsl_vector_get(eval,0), gsl_vector_get(eval,1), gsl_vector_get(eval,2);

    // clean up
    gsl_eigen_symm_free (w);
    gsl_vector_free     (eval);
}

/** Eigenvalues and eigenvectors. NOTE: This function changes the matrix M. */
inline void Eig (Mat3_t & M, Vec3_t & L, Vec3_t & V0, Vec3_t & V1, Vec3_t & V2, bool SortAsc=false, bool SortDesc=false)
{
    // calculate
    gsl_matrix_view m = gsl_matrix_view_array (M.data(), 3, 3);
    gsl_vector * eval = gsl_vector_alloc      (3);
    gsl_matrix * evec = gsl_matrix_alloc      (3, 3);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
    gsl_eigen_symmv (&m.matrix, eval, evec, w);

    // sort
    if (SortAsc)  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    if (SortDesc) gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    // eigenvalues
    L = gsl_vector_get(eval,0), gsl_vector_get(eval,1), gsl_vector_get(eval,2);

    // eigenvectors
    gsl_vector_view ev = gsl_matrix_column (evec,0);
    V0 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);
    ev = gsl_matrix_column (evec,1);
    V1 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);
    ev = gsl_matrix_column (evec,2);
    V2 = gsl_vector_get    (&ev.vector,0), gsl_vector_get(&ev.vector,1), gsl_vector_get(&ev.vector,2);

    // clean up
    gsl_eigen_symmv_free (w);
    gsl_vector_free      (eval);
    gsl_matrix_free      (evec);
}

/** Norm. */
inline double Norm (Vec3_t const & V)
{
    return sqrt(blitz::dot(V,V));
}

/** Dyadic product. */
inline void Dyad (Vec3_t const & A, Vec3_t const & B, Mat3_t & M)
{
    M(0,0)=A(0)*B(0);  M(0,1)=A(0)*B(1);  M(0,2)=A(0)*B(2);
    M(1,0)=A(1)*B(0);  M(1,1)=A(1)*B(1);  M(1,2)=A(1)*B(2);
    M(2,0)=A(2)*B(0);  M(2,1)=A(2)*B(1);  M(2,2)=A(2)*B(2);
}

/** Dyadic product multiplied by s. */
inline void Dyad (double s, Vec3_t const & A, Vec3_t const & B, Mat3_t & M)
{
    M(0,0)=s*A(0)*B(0);  M(0,1)=s*A(0)*B(1);  M(0,2)=s*A(0)*B(2);
    M(1,0)=s*A(1)*B(0);  M(1,1)=s*A(1)*B(1);  M(1,2)=s*A(1)*B(2);
    M(2,0)=s*A(2)*B(0);  M(2,1)=s*A(2)*B(1);  M(2,2)=s*A(2)*B(2);
}

/** Left multiplication. {B} = {A}^T*[M].  */
inline void Mult (Vec3_t const & A, Mat3_t const & M, Vec3_t & B)
{
    B(0) = A(0)*M(0,0) + A(1)*M(1,0) + A(2)*M(2,0);
    B(1) = A(0)*M(0,1) + A(1)*M(1,1) + A(2)*M(2,1);
    B(2) = A(0)*M(0,2) + A(1)*M(1,2) + A(2)*M(2,2);
}

/** Right multiplication. {B} = [M]*{A}.  */
inline void Mult (Mat3_t const & M, Vec3_t const & A, Vec3_t & B)
{
    B(0) = A(0)*M(0,0) + A(1)*M(0,1) + A(2)*M(0,2);
    B(1) = A(0)*M(1,0) + A(1)*M(1,1) + A(2)*M(1,2);
    B(2) = A(0)*M(2,0) + A(1)*M(2,1) + A(2)*M(2,2);
}

/** Matrix multiplication. */
inline void Mult (Mat3_t const & A, Mat3_t const & B, Mat3_t & M)
{
    M(0,0)=A(0,2)*B(2,0)+A(0,1)*B(1,0)+A(0,0)*B(0,0);  M(0,1)=A(0,2)*B(2,1)+A(0,1)*B(1,1)+A(0,0)*B(0,1);  M(0,2)=A(0,2)*B(2,2)+A(0,1)*B(1,2)+A(0,0)*B(0,2);
    M(1,0)=A(1,2)*B(2,0)+B(1,0)*A(1,1)+B(0,0)*A(1,0);  M(1,1)=A(1,2)*B(2,1)+A(1,1)*B(1,1)+B(0,1)*A(1,0);  M(1,2)=A(1,2)*B(2,2)+A(1,1)*B(1,2)+B(0,2)*A(1,0);
    M(2,0)=A(2,2)*B(2,0)+B(1,0)*A(2,1)+B(0,0)*A(2,0);  M(2,1)=B(2,1)*A(2,2)+B(1,1)*A(2,1)+B(0,1)*A(2,0);  M(2,2)=A(2,2)*B(2,2)+B(1,2)*A(2,1)+B(0,2)*A(2,0);
}

/** Clear vector. */
inline void set_to_zero (Vec3_t & V)
{
    V = 0.0, 0.0, 0.0;
}

/** Clear matrix. */
inline void set_to_zero (Mat3_t & M)
{
    M = 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0;
}

// Constants
namespace OrthoSys
{
    Vec3_t O;        ///< Origin
    Vec3_t e0,e1,e2; ///< Basis
    Mat3_t I;        ///< Identity

    int __init_ortho_sys()
    {
        O  = 0.0, 0.0, 0.0;
        e0 = 1.0, 0.0, 0.0;
        e1 = 0.0, 1.0, 0.0;
        e2 = 0.0, 0.0, 1.0;
        I  = 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0;
        return 0.0;
    }

    int __dummy_init_ortho_sys = __init_ortho_sys();
}

#ifdef USE_BOOST_PYTHON
inline Vec3_t Tup2Vec3 (BPy::tuple const & T3)
{
    return Vec3_t(BPy::extract<double>(T3[0])(), BPy::extract<double>(T3[1])(), BPy::extract<double>(T3[2])());
}
#endif


//////////////////////////////////////////////////////////////////////////////////////////// Conversions /////////


/** Creates the matrix representation of 2nd order symmetric tensor Ten (Mandel's representation). */
inline void Ten2Mat (Vec_t const & Ten, Mat3_t & Mat)
{
    // matrix of tensor
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,     0.0,
              Ten(3)/Util::SQ2,  Ten(1),               0.0,
                           0.0,               0.0,  Ten(2);
    }
    else if (ncp==6)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,  Ten(5)/Util::SQ2,
              Ten(3)/Util::SQ2,  Ten(1),            Ten(4)/Util::SQ2,
              Ten(5)/Util::SQ2,  Ten(4)/Util::SQ2,  Ten(2);
        
    }
    else throw new Fatal("matvec.h::Ten2Mat: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** (TinyMatrix) Creates the matrix representation of 2nd order symmetric tensor Ten (Mandel's representation). */
inline void Ten2Mat (Vec_t const & Ten, Mat_t & Mat)
{
    // matrix of tensor
    Mat.change_dim (3,3);
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,     0.0,
              Ten(3)/Util::SQ2,  Ten(1),               0.0,
                           0.0,               0.0,  Ten(2);
    }
    else if (ncp==6)
    {
        Mat = Ten(0),            Ten(3)/Util::SQ2,  Ten(5)/Util::SQ2,
              Ten(3)/Util::SQ2,  Ten(1),            Ten(4)/Util::SQ2,
              Ten(5)/Util::SQ2,  Ten(4)/Util::SQ2,  Ten(2);
        
    }
    else throw new Fatal("matvec.h::Ten2Mat: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Creates a 2nd order symmetric tensor Ten (using Mandel's representation) based on its matrix components. */
inline void Mat2Ten (Mat3_t const & Mat, Vec_t & Ten, size_t NumComponents, bool CheckSymmetry=true, double Tol=1.0e-10)
{
    if (CheckSymmetry)
    {
        if (fabs(Mat(0,1)-Mat(1,0))>Tol) throw new Fatal("matvec.h::Mat2Ten: Matrix is not symmetric");
        if (fabs(Mat(1,2)-Mat(2,1))>Tol) throw new Fatal("matvec.h::Mat2Ten: Matrix is not symmetric");
        if (fabs(Mat(2,0)-Mat(0,2))>Tol) throw new Fatal("matvec.h::Mat2Ten: Matrix is not symmetric");
    }

    // tensor from matrix
    Ten.change_dim (NumComponents);
    if (NumComponents==4)
    {
        Ten = Mat(0,0), Mat(1,1), Mat(2,2), Util::SQ2*Mat(0,1);
    }
    else if (NumComponents==6)
    {
        Ten = Mat(0,0), Mat(1,1), Mat(2,2), Util::SQ2*Mat(0,1), Util::SQ2*Mat(1,2), Util::SQ2*Mat(2,0);
    }
    else throw new Fatal("matvec.h::Mat2Ten: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

#ifdef HAS_TENSORS

/** Converts a Vector to 1st order Tensor. */
inline void Vec2Tensor (Vec_t const & V, Ten1_t & Vec)
{
    if (size(V)!=3) throw new Fatal("matvec.h::Ten1Tensor: Vector (Vec_t) must have size equal to 3");
    Vec = V(0), V(1), V(2);
}

/** Converts a TinyVector to 1st order Tensor. */
inline void Vec2Tensor (Vec3_t const & V, Ten1_t & Vec)
{
    Vec = V(0), V(1), V(2);
}

/** Converts a 1st order Tensor to Vector. */
inline void Tensor2Vec (Ten1_t const & V, Vec_t & Vec)
{
    Vec.change_dim (3);
    Vec = V[0], V[1], V[2];
}

/** Converts a 1st order Tensor to TinyVector. */
inline void Tensor2Vec (Ten1_t const & V, Vec3_t & Vec)
{
    Vec = V[0], V[1], V[2];
}

/** Converts a 2nd order Tensor to TinyMatrix. */
inline void Tensor2Mat (Ten2_t const & T, Mat3_t & M)
{
    M = T[0][0],  T[0][1],  T[0][2],
        T[1][0],  T[1][1],  T[1][2],
        T[2][0],  T[2][1],  T[2][2];
}

/** Converts a TinyMatrix to 2nd order Tensor. */
inline void Mat2Tensor (Mat3_t const & M, Ten2_t & T)
{
    T = M(0,0),  M(0,1),  M(0,2),
        M(1,0),  M(1,1),  M(1,2),
        M(2,0),  M(2,1),  M(2,2);
}

/** Converts a 2nd order symmetric tensor Ten (using Mandel's basis) to 2nd order Tensor. */
inline void Ten2Tensor (Vec_t const & Ten, Ten2_t & T)
{
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        T = Ten(0),            Ten(3)/Util::SQ2,     0.0,
            Ten(3)/Util::SQ2,  Ten(1),               0.0,
                         0.0,               0.0,  Ten(2);
    }
    else if (ncp==6)
    {
        T = Ten(0),            Ten(3)/Util::SQ2,  Ten(5)/Util::SQ2,
            Ten(3)/Util::SQ2,  Ten(1),            Ten(4)/Util::SQ2,
            Ten(5)/Util::SQ2,  Ten(4)/Util::SQ2,  Ten(2);
    }
    else throw new Fatal("matvec.h::Ten2Tensor: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation. NumComponents==%d is invalid",ncp);
}

/** Converts a 2nd order Tensor to tensor in Mandel's basis. */
inline void Tensor2Ten (Ten2_t const & T, Vec_t & Ten, size_t NumComponents, bool CheckSymmetry=true, double Tol=1.0e-10)
{
    if (CheckSymmetry)
    {
        bool error = false;
        if (fabs(T[0][1]-T[1][0])>Tol) error = true;
        if (fabs(T[1][2]-T[2][1])>Tol) error = true;
        if (fabs(T[2][0]-T[0][2])>Tol) error = true;
        if (error)
        {
            Mat3_t mT;
            Tensor2Mat (T, mT);
            std::ostringstream oss;
            oss << "T = \n" << PrintMatrix(mT);
            throw new Fatal("matvec.h::Tensor2Ten: Tensor is not symmetric\n%s",oss.str().c_str());
        }
    }

    Ten.change_dim (NumComponents);
    if (NumComponents==4)
    {
        Ten = T[0][0], T[1][1], T[2][2], T[0][1]*Util::SQ2;
    }
    else if (NumComponents==6)
    {
        Ten = T[0][0], T[1][1], T[2][2], T[0][1]*Util::SQ2, T[1][2]*Util::SQ2, T[2][0]*Util::SQ2;
    }
    else throw new Fatal("matvec.h::Tensor2Ten: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation. NumComponents==%d is invalid",NumComponents);
}

/** Converts a 4th order Tensor to tensor in Mandel's basis. */
inline void Tensor2Ten (Ten4_t const & T, Mat_t & Ten, size_t NumComponents, bool CheckSymmetry=true, double Tol=1.0e-10)
{
    if (!(NumComponents==4 || NumComponents==6)) throw new Fatal("matvec.h::Tensor2Ten: This method is only available for 4th order symmetric tensors with either 4x4 or 6x6 components according to Mandel's representation. NumComponents==%d is invalid",NumComponents);
    Ten.change_dim (NumComponents,NumComponents);
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
    for (size_t k=0; k<3; ++k)
    for (size_t l=0; l<3; ++l)
    {
        if (CheckSymmetry)
        {
            if (fabs(T[i][j][k][l] - T[i][j][l][k])>Tol) throw new Fatal("matvec.h::Tensor2Ten: 4ht order Tensor is not super symmetric: ijkl != ijlk");
            if (fabs(T[i][j][k][l] - T[j][i][k][l])>Tol) throw new Fatal("matvec.h::Tensor2Ten: 4ht order Tensor is not super symmetric: ijkl != jikl");
            if (fabs(T[i][j][k][l] - T[k][l][i][j])>Tol) throw new Fatal("matvec.h::Tensor2Ten: 4ht order Tensor is not super symmetric: ijkl != klij");
        }
        size_t I = TensorsLib::Tensor4ToMandel_i[i][j][k][l];
        size_t J = TensorsLib::Tensor4ToMandel_j[i][j][k][l];
        double a = (k==l ? 1.0 : Util::SQ2);
        double b = (i==j ? 1.0 : Util::SQ2);
        if (I<NumComponents && J<NumComponents) Ten(I,J) = T[i][j][k][l]*a*b;
    }
}

#endif


//////////////////////////////////////////////////////////////////////////////////////// Derivatives /////////////


/** Derivative of unit vector. */
inline void UnitVecDeriv (Vec_t const & n, Vec_t & nu, Mat_t & dnudn, double Tol=1.0e-8)
{
    double norm_n = Norm(n);
    Mat_t I;
    Identity (size(n), I);
    if (norm_n>Tol)
    {
        nu = n / norm_n;
        Mat_t nu_dy_nu;
        Dyad (nu, nu, nu_dy_nu);
        dnudn = I/norm_n - nu_dy_nu/norm_n;
    }
    else
    {
        nu    = 1./Util::SQ3, 1./Util::SQ3, 1./Util::SQ3;
        dnudn = I;
    }
}

/** Derivative of unit vector. */
inline void UnitVecDeriv (Vec3_t const & n, Vec3_t & nu, Mat3_t & dnudn, double Tol=1.0e-8)
{
    double norm_n = Norm(n);
    if (norm_n>Tol)
    {
        double s = 1.0/norm_n;
        nu = n / norm_n;
        dnudn(0,0)=s-s*nu(0)*nu(0);   dnudn(0,1)= -s*nu(0)*nu(1);   dnudn(0,2)= -s*nu(0)*nu(2);
        dnudn(1,0)= -s*nu(1)*nu(0);   dnudn(1,1)=s-s*nu(1)*nu(1);   dnudn(1,2)= -s*nu(1)*nu(2);
        dnudn(2,0)= -s*nu(2)*nu(0);   dnudn(2,1)= -s*nu(2)*nu(1);   dnudn(2,2)=s-s*nu(2)*nu(2);
    }
    else
    {
        nu    = 1./Util::SQ3, 1./Util::SQ3, 1./Util::SQ3;
        dnudn = 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
    }
}

#ifdef HAS_TENSORS

/** Derivative of unit vector (saves results in Tensor). */
inline void UnitVecDeriv (Vec3_t const & n, Vec3_t & nu, Ten2_t & dnudn, double Tol=1.0e-8)
{
    double norm_n = Norm(n);
    if (norm_n>Tol)
    {
        double s = 1.0/norm_n;
        nu = n / norm_n;
        dnudn[0][0]=s-s*nu(0)*nu(0);   dnudn[0][1]= -s*nu(0)*nu(1);   dnudn[0][2]= -s*nu(0)*nu(2);
        dnudn[1][0]= -s*nu(1)*nu(0);   dnudn[1][1]=s-s*nu(1)*nu(1);   dnudn[1][2]= -s*nu(1)*nu(2);
        dnudn[2][0]= -s*nu(2)*nu(0);   dnudn[2][1]= -s*nu(2)*nu(1);   dnudn[2][2]=s-s*nu(2)*nu(2);
    }
    else
    {
        nu    = 1./Util::SQ3, 1./Util::SQ3, 1./Util::SQ3;
        dnudn = 1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
    }
}

#endif


///////////////////////////////////////////////////////////////////////////////////////////// Tensors ////////////


// Assymmetric 3D tensor
typedef blitz::TinyVector<double,9> ATensor2; ///< Assymmetric 2nd order tensor: for deformation gradient (F), velocity gradient (l), ...

/** Symmetric part of a tensor multiplied by a. S = a*sym(T) = a * 0.5*(T + trn(T)). */
inline void Sym (double a, ATensor2 const & T,  size_t NCp, Vec_t & S)
{
    S.change_dim (NCp);
    if (NCp==4)
    {
        S(0) = a*T(0);
        S(1) = a*T(1);
        S(2) = a*T(2);
        S(3) = a*0.5*(T(3)+T(6))*Util::SQ2;
    }
    else if (NCp==6)
    {
        S(0) = a*T(0);
        S(1) = a*T(1);
        S(2) = a*T(2);
        S(3) = a*0.5*(T(3)+T(6))*Util::SQ2;
        S(4) = a*0.5*(T(4)+T(7))*Util::SQ2;
        S(5) = a*0.5*(T(5)+T(8))*Util::SQ2;
    }
    else throw new Fatal("matvec.h::Sym: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation. NCp=%zd is invalid",NCp);
}

/** Dot product: C=A*B  =>  C(i,j)=A(i,k)*B(k,j) */
inline void Dot (ATensor2 const & A, ATensor2 const & B,  ATensor2 & C)
{
    C(0) = A(0)*B(0)+A(3)*B(6)+A(5)*B(8);
    C(3) = A(0)*B(3)+A(3)*B(1)+A(5)*B(7);
    C(5) = A(0)*B(5)+A(3)*B(4)+A(5)*B(2);
    C(6) = A(6)*B(0)+A(1)*B(6)+A(4)*B(8);
    C(1) = A(6)*B(3)+A(1)*B(1)+A(4)*B(7);
    C(4) = A(6)*B(5)+A(1)*B(4)+A(4)*B(2);
    C(8) = A(8)*B(0)+A(7)*B(6)+A(2)*B(8);
    C(7) = A(8)*B(3)+A(7)*B(1)+A(2)*B(7);
    C(2) = A(8)*B(5)+A(7)*B(4)+A(2)*B(2);
}

/** Add R to operation involving skew and symm tensors: R += w*S - S*w, in which w = skew(L) = 0.5*(L-L^T). */
inline void AddSkewTimesOp1 (ATensor2 const & L, Vec_t const & S,  Vec_t & R, double Tol=1.0e-14)
{
    size_t ncp = size(R);
    if (ncp==4)
    {
        R(0) += -(Util::SQ2*S(3)*L(6)-Util::SQ2*L(3)*S(3))/2.0;
        R(1) +=  (Util::SQ2*S(3)*L(6)-Util::SQ2*L(3)*S(3))/2.0;
        R(2) +=  0.0;
        R(3) += -((S(1)-S(0))*L(6)+(S(0)-S(1))*L(3))/Util::SQ2;
        double R4 = (Util::SQ2*S(3)*L(8)+(2*S(1)-2*S(2))*L(7)-Util::SQ2*S(3)*L(5)+(2*S(2)-2*S(1))*L(4))/(2.0*Util::SQ2);
        double R5 = -((2*S(2)-2*S(0))*L(8)-Util::SQ2*S(3)*L(7)+(2*S(0)-2*S(2))*L(5)+Util::SQ2*S(3)*L(4))/(2.0*Util::SQ2);
        if (fabs(R4)>Tol) throw new Fatal("matvec.h::SkewTimesOp1: R4=%g is not zero => Tensor cannot be represented by 4 components only",R4);
        if (fabs(R5)>Tol) throw new Fatal("matvec.h::SkewTimesOp1: R5=%g is not zero => Tensor cannot be represented by 5 components only",R5);
    }
    else if (ncp==6)
    {
        R(0) += -(Util::SQ2*S(5)*L(8)+Util::SQ2*S(3)*L(6)-Util::SQ2*L(5)*S(5)-Util::SQ2*L(3)*S(3))/2.0;
        R(1) += -(Util::SQ2*S(4)*L(7)-Util::SQ2*S(3)*L(6)-Util::SQ2*L(4)*S(4)+Util::SQ2*L(3)*S(3))/2.0;
        R(2) +=  (Util::SQ2*S(5)*L(8)+Util::SQ2*S(4)*L(7)-Util::SQ2*L(5)*S(5)-Util::SQ2*L(4)*S(4))/2.0;
        R(3) += -(Util::SQ2*S(4)*L(8)+Util::SQ2*S(5)*L(7)+(2*S(1)-2*S(0))*L(6)-Util::SQ2*L(4)*S(5)-Util::SQ2*S(4)*L(5)+(2*S(0)-2*S(1))*L(3))/(2.0*Util::SQ2);
        R(4) +=  (Util::SQ2*S(3)*L(8)+(2*S(1)-2*S(2))*L(7)+Util::SQ2*S(5)*L(6)-Util::SQ2*L(3)*S(5)-Util::SQ2*S(3)*L(5)+(2*S(2)-2*S(1))*L(4))/(2.0*Util::SQ2);
        R(5) += -((2*S(2)-2*S(0))*L(8)-Util::SQ2*S(3)*L(7)+Util::SQ2*S(4)*L(6)+(2*S(0)-2*S(2))*L(5)-Util::SQ2*L(3)*S(4)+Util::SQ2*S(3)*L(4))/(2.0*Util::SQ2);
    }
    else throw new Fatal("matvec.h::SkewTimesOp1: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation. NCp=%zd is invalid",ncp);
}

/** Determinant of F. */
inline double Det (ATensor2 const & F)
{
    return -F(3)*(F(2)*F(6)-F(4)*F(8)) + F(5)*(F(6)*F(7)-F(1)*F(8)) + F(0)*(F(1)*F(2)-F(4)*F(7));
}

/** Calculate left Cauchy-Green tensor: b = F*F^T. */
inline void CalcLCauchyGreen (ATensor2 const & F,  size_t NCp, Vec_t & b, double Tol=1.0e-14)
{
    b.change_dim (NCp);
    if (NCp==4)
    {
        b(0) = pow(F(5),2.0)+pow(F(3),2.0)+pow(F(0),2.0);
        b(1) = pow(F(6),2.0)+pow(F(4),2.0)+pow(F(1),2.0);
        b(2) = pow(F(8),2.0)+pow(F(7),2.0)+pow(F(2),2.0);
        b(3) = Util::SQ2*F(0)*F(6)+Util::SQ2*F(4)*F(5)+Util::SQ2*F(1)*F(3);
        double b4 = Util::SQ2*F(6)*F(8)+Util::SQ2*F(1)*F(7)+Util::SQ2*F(2)*F(4);
        double b5 = Util::SQ2*F(0)*F(8)+Util::SQ2*F(3)*F(7)+Util::SQ2*F(2)*F(5);
        if (fabs(b4)>Tol) throw new Fatal("matvec.h::CalcLCauchyGreen: b4=%g is not zero => Tensor cannot be represented by 4 components only",b4);
        if (fabs(b5)>Tol) throw new Fatal("matvec.h::CalcLCauchyGreen: b5=%g is not zero => Tensor cannot be represented by 5 components only",b5);
    }
    else if (NCp==6)
    {
        b(0) = pow(F(5),2.0)+pow(F(3),2.0)+pow(F(0),2.0);
        b(1) = pow(F(6),2.0)+pow(F(4),2.0)+pow(F(1),2.0);
        b(2) = pow(F(8),2.0)+pow(F(7),2.0)+pow(F(2),2.0);
        b(3) = Util::SQ2*F(0)*F(6)+Util::SQ2*F(4)*F(5)+Util::SQ2*F(1)*F(3);
        b(4) = Util::SQ2*F(6)*F(8)+Util::SQ2*F(1)*F(7)+Util::SQ2*F(2)*F(4);
        b(5) = Util::SQ2*F(0)*F(8)+Util::SQ2*F(3)*F(7)+Util::SQ2*F(2)*F(5);
    }
    else throw new Fatal("matvec.h::CalcLCauchyGreen: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation. NCp=%zd is invalid",NCp);
}


/** Multiplication of tensor's matrix by a vector (similar to Cauchy's rule). {t} = [Sig]*{n}. */
inline void Mult (Vec_t const & Sig, Vec3_t const & n, Vec3_t & t)
{
    size_t ncp = size(Sig);
    if (ncp==4)
    {
        t(0) = Sig(0)*n(0)           + Sig(3)*n(1)/Util::SQ2;
        t(1) = Sig(3)*n(0)/Util::SQ2 + Sig(1)*n(1);
        t(2) =                                               + Sig(2)*n(2);
    }
    else if (ncp==6)
    {
        t(0) = Sig(0)*n(0)           + Sig(3)*n(1)/Util::SQ2 + Sig(5)*n(2)/Util::SQ2;
        t(1) = Sig(3)*n(0)/Util::SQ2 + Sig(1)*n(1)           + Sig(4)*n(2)/Util::SQ2;
        t(2) = Sig(5)*n(0)/Util::SQ2 + Sig(4)*n(1)/Util::SQ2 + Sig(2)*n(2);
    }
    else throw new Fatal("matvec.h::Mult: (Cauchy) This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Deviator of 2nd order symmetric tensor Ten. */
inline void Dev (Vec_t const & Ten, Vec_t & DevTen)
{
    double coef = (Ten(0)+Ten(1)+Ten(2))/3.0;
    DevTen     = Ten;
    DevTen(0) -= coef;
    DevTen(1) -= coef;
    DevTen(2) -= coef;
}

/** Trace of 2nd order symmetric tensor Ten. */
inline double Tra (Vec_t const & Ten)
{
    return Ten(0)+Ten(1)+Ten(2);
}

/** 2nd order symmetric tensor Ten raised to the power of 2. */
inline void Pow2 (Vec_t const & Ten, Vec_t & Ten2)
{
    size_t ncp = size(Ten);
    Ten2.change_dim (ncp);
    if (ncp==4)
    {
        Ten2 = Ten(0)*Ten(0) + Ten(3)*Ten(3)/2.0,
               Ten(1)*Ten(1) + Ten(3)*Ten(3)/2.0,
               Ten(2)*Ten(2),
               Ten(0)*Ten(3) + Ten(1)*Ten(3);
    }
    else if (ncp==6)
    {
        Ten2 = Ten(0)*Ten(0)           + Ten(3)*Ten(3)/2.0       + Ten(5)*Ten(5)/2.0,
               Ten(3)*Ten(3)/2.0       + Ten(1)*Ten(1)           + Ten(4)*Ten(4)/2.0,
               Ten(5)*Ten(5)/2.0       + Ten(4)*Ten(4)/2.0       + Ten(2)*Ten(2),
               Ten(0)*Ten(3)           + Ten(3)*Ten(1)           + Ten(5)*Ten(4)/Util::SQ2,
               Ten(3)*Ten(5)/Util::SQ2 + Ten(1)*Ten(4)           + Ten(4)*Ten(2),
               Ten(0)*Ten(5)           + Ten(3)*Ten(4)/Util::SQ2 + Ten(5)*Ten(2);
    }
    else throw new Fatal("matvec.h::Pow2: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Determinant of 2nd order symmetric tensor Ten. */
inline double Det (Vec_t const & Ten)
{
    size_t ncp = size(Ten);
    if (ncp==4)
    {
        return   Ten(0)*Ten(1)*Ten(2)
               - Ten(2)*Ten(3)*Ten(3)/2.0;
    }
    else if (ncp==6)
    {
        return   Ten(0)*Ten(1)*Ten(2) 
               + Ten(3)*Ten(4)*Ten(5)/Util::SQ2
               - Ten(0)*Ten(4)*Ten(4)/2.0
               - Ten(1)*Ten(5)*Ten(5)/2.0
               - Ten(2)*Ten(3)*Ten(3)/2.0;
    }
    else throw new Fatal("matvec.h::Det: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Inverse of 2nd order symmetric tensor T. */
inline void Inv (Vec_t const & T, Vec_t & Ti, double Tol=1.0e-10)
{
    size_t ncp = size(T);
    Ti.change_dim (ncp);
    if (ncp==4)
    {
        double det =  T(0)*T(1)*T(2)
                    - T(2)*T(3)*T(3)/2.0;

        if (fabs(det)<Tol)
        {
            std::ostringstream oss;  oss<<PrintVector(T);
            throw new Fatal("matvec.h::Inv: inverse of 2nd order symmetric tensor failed with null (%g) determinat.\n  T =%s",Tol,oss.str().c_str());
        }

        Ti(0) =  T(1)*T(2)/det;
        Ti(1) =  T(0)*T(2)/det;
        Ti(2) = (T(0)*T(1) - T(3)*T(3)/2.0)/det;
        Ti(3) = -T(2)*T(3)/det;

    }
    else if (ncp==6)
    {
        double det =  T(0)*T(1)*T(2) 
                    + T(3)*T(4)*T(5)/Util::SQ2
                    - T(0)*T(4)*T(4)/2.0
                    - T(1)*T(5)*T(5)/2.0
                    - T(2)*T(3)*T(3)/2.0;

        if (fabs(det)<Tol)
        {
            std::ostringstream oss;  oss<<PrintVector(T);
            throw new Fatal("matvec.h::Inv: inverse of 2nd order symmetric tensor failed with null (%g) determinat.\n  T =%s",Tol,oss.str().c_str());
        }

        Ti(0) = (T(1)*T(2) - T(4)*T(4)/2.0)/det;
        Ti(1) = (T(0)*T(2) - T(5)*T(5)/2.0)/det;
        Ti(2) = (T(0)*T(1) - T(3)*T(3)/2.0)/det;
        Ti(3) = (T(4)*T(5)/Util::SQ2 - T(2)*T(3))/det;
        Ti(4) = (T(3)*T(5)/Util::SQ2 - T(0)*T(4))/det;
        Ti(5) = (T(3)*T(4)/Util::SQ2 - T(1)*T(5))/det;
    }
    else throw new Fatal("matvec.h::Inv: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Characteristic invariants of 2nd order symmetric tensor Ten. */
inline void CharInvs (Vec_t const & Ten, double & I1, double & I2, double & I3)
{
    I1 = Ten(0) + Ten(1) + Ten(2);
    I2 = Ten(0)*Ten(1) + Ten(1)*Ten(2) + Ten(2)*Ten(0) - Ten(3)*Ten(3)/2.0;
    I3 = Ten(0)*Ten(1)*Ten(2) - Ten(2)*Ten(3)*Ten(3)/2.0;
    size_t ncp = size(Ten);
    if (ncp>4)
    {
        I2 += (-Ten(4)*Ten(4)/2.0 - Ten(5)*Ten(5)/2.0);
        I3 += (Ten(3)*Ten(4)*Ten(5)/Util::SQ2 - Ten(0)*Ten(4)*Ten(4)/2.0 - Ten(1)*Ten(5)*Ten(5)/2.0);
    }
}

/** Characteristic invariants of 2nd order symmetric tensor Ten and their derivatives. */
inline void CharInvs (Vec_t const & Ten, double & I1, double & I2, double & I3, Vec_t & dI1dTen, Vec_t & dI2dTen, Vec_t & dI3dTen)
{
    I1 = Ten(0) + Ten(1) + Ten(2);
    I2 = Ten(0)*Ten(1) + Ten(1)*Ten(2) + Ten(2)*Ten(0) - Ten(3)*Ten(3)/2.0;
    I3 = Ten(0)*Ten(1)*Ten(2) - Ten(2)*Ten(3)*Ten(3)/2.0;
    size_t ncp = size(Ten);
    dI1dTen.change_dim (ncp);
    dI2dTen.change_dim (ncp);
    dI3dTen.change_dim (ncp);
    if (ncp>4)
    {
        I2 += (-Ten(4)*Ten(4)/2.0 - Ten(5)*Ten(5)/2.0);
        I3 += (Ten(3)*Ten(4)*Ten(5)/Util::SQ2 - Ten(0)*Ten(4)*Ten(4)/2.0 - Ten(1)*Ten(5)*Ten(5)/2.0);
        dI1dTen = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    }
    else dI1dTen = 1.0, 1.0, 1.0, 0.0;
    Vec_t Ten2(ncp);
    Pow2 (Ten, Ten2);
    dI2dTen = I1*dI1dTen - Ten;
    dI3dTen = Ten2 - I1*Ten + I2*dI1dTen;
}

/** Derivative of the Inverse of A w.r.t A. */
inline void DerivInv (Vec_t const & A, Vec_t & Ai, Mat_t & dInvA_dA, double Tol=1.0e-10)
{
    Inv (A, Ai, Tol);
    size_t ncp = size(A);
    dInvA_dA.change_dim (ncp,ncp);
    if (ncp==4) throw new Fatal("matvec.h: DerivInv: This method doesn't work for ncp==4");
    else if (ncp==6)
    {
        double s = Util::SQ2;
        dInvA_dA =
            -Ai(0)*Ai(0)     , -Ai(3)*Ai(3)/2.  , -Ai(5)*Ai(5)/2.  , -Ai(0)*Ai(3)                      , -(Ai(3)*Ai(5))/s                  , -Ai(0)*Ai(5)                      , 
            -Ai(3)*Ai(3)/2.  , -Ai(1)*Ai(1)     , -Ai(4)*Ai(4)/2.  , -Ai(1)*Ai(3)                      , -Ai(1)*Ai(4)                      , -(Ai(3)*Ai(4))/s                  , 
            -Ai(5)*Ai(5)/2.  , -Ai(4)*Ai(4)/2.  , -Ai(2)*Ai(2)     , -(Ai(4)*Ai(5))/s                  , -Ai(2)*Ai(4)                      , -Ai(2)*Ai(5)                      , 
            -Ai(0)*Ai(3)     , -Ai(1)*Ai(3)     , -(Ai(4)*Ai(5))/s , -Ai(3)*Ai(3)/2.-Ai(0)*Ai(1)       , -(Ai(1)*Ai(5))/s-(Ai(3)*Ai(4))/2. , -(Ai(3)*Ai(5))/2.-(Ai(0)*Ai(4))/s , 
            -(Ai(3)*Ai(5))/s , -Ai(1)*Ai(4)     , -Ai(2)*Ai(4)     , -(Ai(1)*Ai(5))/s-(Ai(3)*Ai(4))/2. , -Ai(4)*Ai(4)/2.-Ai(1)*Ai(2)       , -(Ai(4)*Ai(5))/2.-(Ai(2)*Ai(3))/s , 
            -Ai(0)*Ai(5)     , -(Ai(3)*Ai(4))/s , -Ai(2)*Ai(5)     , -(Ai(3)*Ai(5))/2.-(Ai(0)*Ai(4))/s , -(Ai(4)*Ai(5))/2.-(Ai(2)*Ai(3))/s , -Ai(5)*Ai(5)/2.-Ai(0)*Ai(2)       ;
    }
    else throw new Fatal("matvec.h::DerivInv: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Eigenvalues and eigenvectors of second order tensor in Mandel's basis. */
inline void Eig (Vec_t const & Ten, Vec3_t & L, Vec3_t & V0, Vec3_t & V1, Vec3_t & V2, bool SortAsc=false, bool SortDesc=false)
{
    Mat3_t M;
    Ten2Mat (Ten, M);
    Eig (M, L, V0, V1, V2, SortAsc, SortDesc);
}

/** Eigenprojectors of 2nd order symmetric tensor Ten. */
inline void EigenProj (Vec_t const & Ten, Vec3_t & L, Vec3_t & v0, Vec3_t & v1, Vec3_t & v2, Vec_t & P0, Vec_t & P1, Vec_t & P2, bool SortAsc=false, bool SortDesc=false)
{
    // matrix of tensor
    Mat3_t ten;
    Ten2Mat (Ten, ten);

    // eigen-values and vectors
    Eig (ten, L, v0, v1, v2, SortAsc, SortDesc);

    // eigen-projectors
    size_t ncp = size(Ten);
    P0.change_dim (ncp);
    P1.change_dim (ncp);
    P2.change_dim (ncp);
    if (ncp==4)
    {
        P0 = v0(0)*v0(0), v0(1)*v0(1), v0(2)*v0(2), v0(0)*v0(1)*Util::SQ2;
        P1 = v1(0)*v1(0), v1(1)*v1(1), v1(2)*v1(2), v1(0)*v1(1)*Util::SQ2;
        P2 = v2(0)*v2(0), v2(1)*v2(1), v2(2)*v2(2), v2(0)*v2(1)*Util::SQ2;
    }
    else
    {
        P0 = v0(0)*v0(0), v0(1)*v0(1), v0(2)*v0(2), v0(0)*v0(1)*Util::SQ2, v0(1)*v0(2)*Util::SQ2, v0(2)*v0(0)*Util::SQ2;
        P1 = v1(0)*v1(0), v1(1)*v1(1), v1(2)*v1(2), v1(0)*v1(1)*Util::SQ2, v1(1)*v1(2)*Util::SQ2, v1(2)*v1(0)*Util::SQ2;
        P2 = v2(0)*v2(0), v2(1)*v2(1), v2(2)*v2(2), v2(0)*v2(1)*Util::SQ2, v2(1)*v2(2)*Util::SQ2, v2(2)*v2(0)*Util::SQ2;
    }
}

/** Eigenprojectors of 2nd order symmetric tensor Ten. */
inline void EigenProjAnalytic (Vec_t const & Ten, Vec3_t & L, Vec_t & P0, Vec_t & P1, Vec_t & P2)
{
    // identity tensor
    size_t ncp = size(Ten);
    Vec_t I(ncp);
    set_to_zero(I);
    I(0)=1.0;  I(1)=1.0;  I(2)=1.0;

    // characteristics invariants
    double I1,I2,I3;
    CharInvs (Ten,I1,I2,I3);

    // inverse tensor
    Vec_t Ti;
    Inv (Ten, Ti);

    // eigen-values/projectors
    P0.change_dim (ncp);
    P1.change_dim (ncp);
    P2.change_dim (ncp);
    Vec_t * P[3] = {&P0, &P1, &P2}; // all three eigen projectors
    double alpha = 2.0*sqrt(I1*I1-3.0*I2);
    double theta = acos((2.0*I1*I1*I1-9.0*I1*I2+27.0*I3)/(2.0*pow(I1*I1-3.0*I2,3.0/2.0)));
    //std::cout << "I1 = " << I1 << ",  I2 = " << I2 << ",  I3 = " << I3 << ",  I1*I1-3*I2 = " << I1*I1-3.0*I2 << std::endl;
    //std::cout << "alpha = " << alpha << ",  theta = " << theta << std::endl;
    for (size_t k=0; k<3; ++k)
    {
        L(k) = (I1+alpha*cos((theta+2.0*Util::PI*(1.0+k))/3.0))/3.0;
        if (fabs(L(k))<1.0e-14) throw new Fatal("matvec.h:EigenProjAnalytic: L(%d)=%g must be non-zero",k,L(k));
        double coef = L(k)/(2.0*L(k)*L(k)-L(k)*I1+I3/L(k));
        (*P[k]) = coef*(Ten+(L(k)-I1)*I+(I3/L(k))*Ti);
    }
}

/** Initialize second order identity tensor. */
inline void Calc_I (size_t NCp, Vec_t & I)
{
    I.change_dim(NCp);
    if      (NCp==4) I = 1.0, 1.0, 1.0, 0.0;
    else if (NCp==6) I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    else throw new Fatal("matvec.h::Calc_I: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order identity tensor. */
inline void Calc_IIsym (size_t NCp, Mat_t & IIsym)
{
    IIsym.change_dim(NCp,NCp);
    if (NCp==4)
    {
        IIsym = 1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0;
    }
    else if (NCp==6)
    {
        IIsym = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    }
    else throw new Fatal("matvec.h::Calc_IIsym: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize IdyI fourth order tensor. */
inline void Calc_IdyI (size_t NCp, Mat_t & IdyI)
{
    IdyI.change_dim(NCp,NCp);
    if (NCp==4)
    {
        IdyI = 1.0, 1.0, 1.0, 0.0,
               1.0, 1.0, 1.0, 0.0,
               1.0, 1.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 0.0;
    }
    else if (NCp==6)
    {
        IdyI = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
    else throw new Fatal("matvec.h::Calc_IdyI: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order symmetric-deviatoric tensor. */
inline void Calc_Psd (size_t NCp, Mat_t & Psd)
{
    Psd.change_dim(NCp,NCp);
    if (NCp==4)
    {
        Psd =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0,
              -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0,
              -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0,
                   0.0,      0.0,      0.0, 1.0;
    }
    else if (NCp==6)
    {
        Psd =  2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
              -1.0/3.0,  2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
              -1.0/3.0, -1.0/3.0,  2.0/3.0, 0.0, 0.0, 0.0,
                   0.0,      0.0,      0.0, 1.0, 0.0, 0.0,
                   0.0,      0.0,      0.0, 0.0, 1.0, 0.0,
                   0.0,      0.0,      0.0, 0.0, 0.0, 1.0;
    }
    else throw new Fatal("matvec.h::Calc_Psd: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Initialize fourth order isotropic tensor. */
inline void Calc_Piso (size_t NCp, Mat_t & Piso)
{
    Piso.change_dim(NCp,NCp);
    if (NCp==4)
    {
        Piso = 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
                   0.0,     0.0,     0.0, 0.0;
    }
    else if (NCp==6)
    {
        Piso = 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
               1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0,
                   0.0,     0.0,     0.0, 0.0, 0.0, 0.0;
    }
    else throw new Fatal("matvec.h::Calc_Piso: This method is only available for 2nd order symmetric tensors with either 4 or 6 components according to Mandel's representation");
}

/** Derivatives of the eigenprojectors w.r.t its defining tensor. */
inline void EigenProjDerivs (Vec_t const & A, Vec3_t & L, Vec3_t & v0,   Vec3_t & v1,   Vec3_t & v2,
                                                          Vec_t &  P0,   Vec_t &  P1,   Vec_t &  P2,
                                                          Mat_t & dP0dA, Mat_t & dP1dA, Mat_t & dP2dA,
                             double Pertubation=1.0e-7, double Tol=1.0e-14, bool SortAsc=false, bool SortDesc=false)
{
    // eigenprojectors
    EigenProj (A, L, v0, v1, v2, P0, P1, P2, SortAsc, SortDesc);

    // check eigenvalues
    if (fabs(L(0))<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));
    if (fabs(L(1))<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));
    if (fabs(L(2))<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));

    // apply pertubation
    bool pert[3] = { false, false, false };
    if (fabs(L(0)-L(1))<Pertubation) { L(0)+=Pertubation;  L(1)-=Pertubation;  pert[0]=true; }
    if (fabs(L(1)-L(2))<Pertubation) { L(1)+=Pertubation;  L(2)-=Pertubation;  pert[1]=true; }
    if (fabs(L(2)-L(0))<Pertubation) { L(2)+=Pertubation;  L(0)-=Pertubation;  pert[2]=true; }

    // characteristic invariants
    double I1 = L(0)+L(1)+L(2);
    double I3 = L(0)*L(1)*L(2);
    double a[3] = {2.*L(0)*L(0)-I1*L(0)+I3/L(0),
                   2.*L(1)*L(1)-I1*L(1)+I3/L(1),
                   2.*L(2)*L(2)-I1*L(2)+I3/L(2)};

    // check alphas
    if (fabs(a[0])<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Alpha0 cannot be zero (Tol=%g).\n  L = [%g, %g, %g]\n  I1,I3 = [%g, %g, %g]\n  Alphas = [%g, %g, %g]",Tol,L(0),L(1),L(2), I1,I3, a[0],a[1],a[2]);
    if (fabs(a[1])<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Alpha1 cannot be zero (Tol=%g).\n  L = [%g, %g, %g]\n  I1,I3 = [%g, %g, %g]\n  Alphas = [%g, %g, %g]",Tol,L(0),L(1),L(2), I1,I3, a[0],a[1],a[2]);
    if (fabs(a[1])<Tol) throw new Fatal("matvec.h::EigenProjDerivs: Alpha2 cannot be zero (Tol=%g).\n  L = [%g, %g, %g]\n  I1,I3 = [%g, %g, %g]\n  Alphas = [%g, %g, %g]",Tol,L(0),L(1),L(2), I1,I3, a[0],a[1],a[2]);

    // inverse tensor and derivative of inverse
    size_t ncp = size(A);
    Vec_t Ai;
    Mat_t dInvAdA, IIsym;
    DerivInv   (A, Ai, dInvAdA, Tol);
    Calc_IIsym (ncp, IIsym);

    // dyadic between projectors
    Array<Mat_t> PdyP(3);
    Dyad (P0,P0, PdyP[0]);
    Dyad (P1,P1, PdyP[1]);
    Dyad (P2,P2, PdyP[2]);

    // constants
    double c[3][3] = {{I3/(a[0]*L(0)*L(0))-L(0)/a[0],  I3/(a[0]*L(1)*L(1))-L(0)/a[0],  I3/(a[0]*L(2)*L(2))-L(0)/a[0]},
                      {I3/(a[1]*L(0)*L(0))-L(1)/a[1],  I3/(a[1]*L(1)*L(1))-L(1)/a[1],  I3/(a[1]*L(2)*L(2))-L(1)/a[1]},
                      {I3/(a[2]*L(0)*L(0))-L(2)/a[2],  I3/(a[2]*L(1)*L(1))-L(2)/a[2],  I3/(a[2]*L(2)*L(2))-L(2)/a[2]}};

    // derivative of projector
    dP0dA.change_dim (ncp,ncp);
    dP1dA.change_dim (ncp,ncp);
    dP2dA.change_dim (ncp,ncp);
    dP0dA = (L(0)/a[0])*IIsym + (I3/a[0])*dInvAdA + c[0][0]*PdyP[0] + c[0][1]*PdyP[1] + c[0][2]*PdyP[2];
    dP1dA = (L(1)/a[1])*IIsym + (I3/a[1])*dInvAdA + c[1][0]*PdyP[0] + c[1][1]*PdyP[1] + c[1][2]*PdyP[2];
    dP2dA = (L(2)/a[2])*IIsym + (I3/a[2])*dInvAdA + c[2][0]*PdyP[0] + c[2][1]*PdyP[1] + c[2][2]*PdyP[2];

    // remove pertubation
    if (pert[0]) { L(0)-=Pertubation;  L(1)+=Pertubation; }
    if (pert[1]) { L(1)-=Pertubation;  L(2)+=Pertubation; }
    if (pert[2]) { L(2)-=Pertubation;  L(0)+=Pertubation; }
}

/** Derivatives of the eigenvectors of Sig w.r.t Sig (a third order tensor).  NOTE: the eigevalues and eigenvectors must be given as input. */
#ifdef HAS_TENSORS
inline void EigenVecDerivs (Vec3_t const &  L, Ten1_t const &  v0,     Ten1_t const &  v1,     Ten1_t const &  v2,
                                               Ten3_t       & dv0dSig, Ten3_t       & dv1dSig, Ten3_t       & dv2dSig,
                            double Pertubation=1.0e-7, double Tol=1.0e-10)
{
    // check eigenvalues
    if (fabs(L(0))<Tol) throw new Fatal("matvec.h::EigenVecDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));
    if (fabs(L(1))<Tol) throw new Fatal("matvec.h::EigenVecDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));
    if (fabs(L(2))<Tol) throw new Fatal("matvec.h::EigenVecDerivs: Principal values cannot be zero (Tol=%g).\n L = [%g, %g, %g]",Tol,L(0),L(1),L(2));

    // apply pertubation
    Vec3_t l(L);
    if (fabs(l(0)-l(1))<Pertubation) { l(0)+=Pertubation;  l(1)-=Pertubation; }
    if (fabs(l(1)-l(2))<Pertubation) { l(1)+=Pertubation;  l(2)-=Pertubation; }
    if (fabs(l(2)-l(0))<Pertubation) { l(2)+=Pertubation;  l(0)-=Pertubation; }

    // derivatives
    double a0 = 0.5/(l(0)-l(1));   double b0 = 0.5/(l(0)-l(2));
    double a1 = 0.5/(l(1)-l(2));   double b1 = 0.5/(l(1)-l(0));
    double a2 = 0.5/(l(2)-l(0));   double b2 = 0.5/(l(2)-l(1));
    dv0dSig = ((v1&v0&v1) + (v1&v1&v0))*a0  +  ((v2&v0&v2) + (v2&v2&v0))*b0;
    dv1dSig = ((v2&v1&v2) + (v2&v2&v1))*a1  +  ((v0&v1&v0) + (v0&v0&v1))*b1;
    dv2dSig = ((v0&v2&v0) + (v0&v0&v2))*a2  +  ((v1&v2&v1) + (v1&v1&v2))*b2;
}
#endif


////////////////////////////////////////////////////////////////////////////////////////// Invariants ////////////


// Cambridge invariants
inline double Calc_pcam  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/3.0; }
inline double Calc_ev    (Vec_t const & Eps) { return   Eps(0)+Eps(1)+Eps(2);      }
inline double Calc_qcam  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(2.0); }
inline double Calc_ed    (Vec_t const & Eps) { double m = (size(Eps)>4 ? pow(Eps(4),2.0)+pow(Eps(5),2.0) : 0.0); return sqrt(pow(Eps(0)-Eps(1),2.0) + pow(Eps(1)-Eps(2),2.0) + pow(Eps(2)-Eps(0),2.0) + 3.0*(pow(Eps(3),2.0)+m))*(sqrt(2.0)/3.0); }

// Octahedral invariants
inline double Calc_poct  (Vec_t const & Sig) { return -(Sig(0)+Sig(1)+Sig(2))/sqrt(3.0); }
inline double Calc_evoct (Vec_t const & Eps) { return  (Eps(0)+Eps(1)+Eps(2))/sqrt(3.0); }
inline double Calc_qoct  (Vec_t const & Sig) { double m = (size(Sig)>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0); return sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(3.0); }
inline double Calc_edoct (Vec_t const & Eps) { double m = (size(Eps)>4 ? pow(Eps(4),2.0)+pow(Eps(5),2.0) : 0.0); return sqrt(pow(Eps(0)-Eps(1),2.0) + pow(Eps(1)-Eps(2),2.0) + pow(Eps(2)-Eps(0),2.0) + 3.0*(pow(Eps(3),2.0)+m))/sqrt(3.0); }

// Octahedral invariants of Sig. */
inline void Calc_pqoct (Vec_t const & Sig, double & p, double & q) { p = Calc_poct(Sig);  q = Calc_qoct(Sig); }


/** Octahedral invariants of Sig. */
inline void OctInvs (Vec_t const & Sig, double & p, double & q, double & t, double qTol=1.0e-8)
{
    size_t ncp = size(Sig);
    Vec_t s;
    Dev (Sig, s);
    double m = (ncp>4 ? pow(Sig(4),2.0)+pow(Sig(5),2.0) : 0.0);
    p = -(Sig(0)+Sig(1)+Sig(2))/Util::SQ3;
    q = sqrt(pow(Sig(0)-Sig(1),2.0) + pow(Sig(1)-Sig(2),2.0) + pow(Sig(2)-Sig(0),2.0) + 3.0*(pow(Sig(3),2.0)+m))/sqrt(3.0);
    t = 0.0;
    if (q>qTol)
    {
        double det_s = Det(s);
        t = -3.0*Util::SQ6*det_s/(q*q*q);
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Octahedral invariants of Sig (and deviator s). */
inline void OctInvs (Vec_t const & Sig, double & p, double & q, Vec_t & s, double qTol=1.0e-8, bool ApplyPertub=false)
{
    q = Calc_qoct (Sig);
    if (q>qTol || !ApplyPertub)
    {
        p = Calc_poct (Sig);
        Dev (Sig, s);
    }
    else
    {
        Vec_t sig(Sig);
        if (sig(3)<0.0) sig(3) -= qTol*Util::SQ2;
        else            sig(3) += qTol*Util::SQ2;
        q = Calc_qoct (sig);
        p = Calc_poct (sig);
        Dev (sig, s);
        if (q<=qTol)
        {
            std::ostringstream oss;
            oss << "Sig = " << PrintVector(Sig);
            oss << "sig = " << PrintVector(sig);
            throw new Fatal("matvec.h::OctInvs: __internal_error__ Pertubation for q<=qTol failed (q=%g, qTol=%g)\n%s",q,qTol,oss.str().c_str());
        }
    }
}

/** Octahedral invariants of Sig (and derivatives). */
inline void OctInvs (Vec_t const & Sig, double & p, double & q, double & t, double & th, Vec_t & s, double qTol=1.0e-8, Vec_t * dthdSig=NULL, bool ApplyPertub=false)
{
    OctInvs (Sig, p, q, s, qTol, ApplyPertub);
    t  = 0.0;
    th = asin(t)/3.0;
    double q3 = pow(q,3.0);
    if (q3>qTol)
    {
        double det_s = Det(s);
        t = -3.0*Util::SQ6*det_s/q3;
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
        th = asin(t)/3.0;
        if (dthdSig!=NULL)
        {
            Vec_t ss, devss;
            Pow2 (s, ss);
            Dev  (ss, devss);
            if (t>-0.999 && t<0.999) (*dthdSig) = (-1.0/(q*q*cos(3.0*th)))*(t*s + (Util::SQ6/q)*devss);
            else set_to_zero ((*dthdSig));
        }
    }
    else if (dthdSig!=NULL) set_to_zero ((*dthdSig));
}

/** Octahedral invariants of Sig (L = principal values). */
inline void OctInvs (Vec3_t const & L, double & p, double & q, double & t, double qTol=1.0e-8)
{
    p = -(L(0)+L(1)+L(2))/Util::SQ3;
    q = sqrt(pow(L(0)-L(1),2.0) + pow(L(1)-L(2),2.0) + pow(L(2)-L(0),2.0))/sqrt(3.0);
    t = 0.0;
    if (q>qTol)
    {
        Vec3_t S((2.*L(0)-L(1)-L(2))/3.,
                 (2.*L(1)-L(2)-L(0))/3.,
                 (2.*L(2)-L(0)-L(1))/3.);
        t = -3.0*Util::SQ6*S(0)*S(1)*S(2)/pow(q,3.0);
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}

/** Octahedral invariants of Sig (L = principal values). With derivatives */
inline void OctInvs (Vec3_t const & L, double & p, double & q, double & t, Vec3_t & dpdL, Vec3_t & dqdL, Vec3_t & dtdL, double qTol=1.0e-8)
{
    Vec3_t one(1.0,1.0,1.0), s;
    p    = -(L(0)+L(1)+L(2))/Util::SQ3;
    q    = sqrt(pow(L(0)-L(1),2.0) + pow(L(1)-L(2),2.0) + pow(L(2)-L(0),2.0))/sqrt(3.0);
    t    = 0.0;
    s    = L - ((L(0)+L(1)+L(2))/3.0)*one;
    dpdL = (-1.0/Util::SQ3)*one;
    dqdL = 0.0, 0.0, 0.0;
    if (q>qTol)
    {
        double q3 = q*q*q;
        double q5 = q3*q*q;
        double l  = (L(0)-L(1))*(L(1)-L(2))*(L(2)-L(0));
        Vec3_t B(L(2)-L(1), L(0)-L(2), L(1)-L(0));
        t    = -3.0*Util::SQ6*s(0)*s(1)*s(2)/q3;
        dqdL = (1.0/q)*s;
        dtdL = (-Util::SQ6*l/q5)*B;
        if (t<=-1.0) t = -1.0;
        if (t>= 1.0) t =  1.0;
    }
}


inline void OctDerivs (Vec3_t const & L, double & p, double & q, double & t, Mat3_t & dpqthdL, double qTol=1.0e-8)
{
    OctInvs (L, p,q,t, qTol);
    if (q>qTol)
    {

        Vec3_t S((2.*L(0)-L(1)-L(2))/3., 
                 (2.*L(1)-L(2)-L(0))/3.,
                 (2.*L(2)-L(0)-L(1))/3.);
        Vec3_t devSS((2.*S(1)*S(2) - S(2)*S(0) - S(0)*S(1))/3.,
                     (2.*S(2)*S(0) - S(1)*S(2) - S(0)*S(1))/3.,
                     (2.*S(0)*S(1) - S(1)*S(2) - S(2)*S(0))/3.);
        double th = asin(t)/3.0;
        double c  = -1.0/(pow(q,3.0)*cos(3.0*th));
        dpqthdL = -1./Util::SQ3,                   -1./Util::SQ3,                   -1./Util::SQ3,
                  S(0)/q,                          S(1)/q,                          S(2)/q,
                  c*Util::SQ6*devSS(0)+c*q*t*S(0), c*Util::SQ6*devSS(1)+c*q*t*S(1), c*Util::SQ6*devSS(2)+c*q*t*S(2);
    }
    else
    {
        dpqthdL = -1./Util::SQ3, -1./Util::SQ3, -1./Util::SQ3,
                   0.,            0.,            0.,
                   0.,            0.,            0.;
    }
}

inline void InvOctDerivs (Vec3_t const & Lsorted, double & p, double & q, double & t, Mat3_t & dLdpqth, double qTol=1.0e-8)
{
    OctInvs (Lsorted, p,q,t, qTol);
    double th = asin(t)/3.0;
    dLdpqth = -1./Util::SQ3, (2.*sin(th-(2.*Util::PI)/3.))/Util::SQ6, (2.*q*cos(th-(2.*Util::PI)/3.))/Util::SQ6,
              -1./Util::SQ3, (2.*sin(th)                 )/Util::SQ6, (2.*q*cos(th)                 )/Util::SQ6,
              -1./Util::SQ3, (2.*sin(th+(2.*Util::PI)/3.))/Util::SQ6, (2.*q*cos(th+(2.*Util::PI)/3.))/Util::SQ6;
}


/////////////////////////////////////////////////////////////////////////////////// Model functions //////////////

/** Calculate K given E and nu. */
inline double Calc_K (double E, double nu) { return E/(3.0*(1.0-2.0*nu)); }

/** Calculate G given E and nu. */
inline double Calc_G (double E, double nu) { return E/(2.0*(1.0+nu)); }

/** Calculate lam given E and nu. */
inline double Calc_lam (double E, double nu) { return E*nu/((1.0+nu)*(1.0-2.0*nu)); }

/** Calculate G given K and nu. */
inline double Calc_G_ (double K, double nu) { return 3.0*(1.0-2.0*nu)*K/(2.0*(1.0+nu)); }

/** Calculate E given K and G. */
inline double Calc_E (double K, double G) { return 9.0*K*G/(3.0*K+G); }

/** Calculate E given K and nu. */
inline double Calc_E_ (double K, double nu) { return 3.0*K*(1.0-2.0*nu); }

/** Calculate nu given K and G. */
inline double Calc_nu (double K, double G) { return (3.0*K-2.0*G)/(6.0*K+2.0*G); }


///////////////////////////////////////////////////////////////////////// Python functions ///////////////////////


#ifdef USE_BOOST_PYTHON

inline void Vec2List (Vec3_t const & V, BPy::list & L)
{
    for (size_t i=0; i<3; ++i) L.append(V(i));
}

inline void Vec2List (Vec_t const & V, BPy::list & L)
{
    for (size_t i=0; i<size(V); ++i) L.append(V(i));
}

inline void List2Vec (BPy::list const & L, Vec_t & V)
{
    size_t m = BPy::len(L);
    if (m>0)
    {
        V.change_dim (m);
        for (size_t i=0; i<m; ++i) V(i) = BPy::extract<double>(L[i])();
    }
}

inline void Mat2List (Mat_t const & M, BPy::list & L)
{
    size_t m = num_rows(M);
    size_t n = num_cols(M);
    for (size_t i=0; i<m; ++i)
    {
        BPy::list row;
        for (size_t j=0; j<n; ++j) row.append(M(i,j));
        L.append(row);
    }
}

inline void List2Mat (BPy::list const & L, Mat_t & M)
{
    size_t m = BPy::len(L);
    if (m>0)
    {
        size_t n = BPy::len(BPy::extract<BPy::list>(L[0])());
        M.change_dim (m,n);
        for (size_t i=0; i<m; ++i)
        {
            BPy::list const & row = BPy::extract<BPy::list>(L[i])();
            for (size_t j=0; j<n; ++j) M(i,j) = BPy::extract<double>(row[j])();
        }
    }
}

inline void PyEigenProjAnalytic (BPy::list const & Ten, BPy::list & L, BPy::list & P0, BPy::list & P1, BPy::list & P2)
{
    Vec3_t l;
    Vec_t  ten, p0, p1, p2;
    List2Vec          (Ten, ten);
    EigenProjAnalytic (ten, l, p0, p1, p2);
    Vec2List          (l, L);
    Vec2List          (p0, P0);
    Vec2List          (p1, P1);
    Vec2List          (p2, P2);
}

#endif

#endif // MECHSYS_MATVEC_H
