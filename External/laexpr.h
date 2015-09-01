/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

/* Expression template evaluation for matrices */

#ifndef MECHSYS_LINALG_LAEXPR_H
#define MECHSYS_LINALG_LAEXPR_H

// STL
#include <iostream>
#include <sstream>

// MechSys
#include "matrix.h"
#include "vector.h"
#include "fatal.h"

extern "C"
{
    // BLAS
    double ddot_(int const *N, double const *X, int const *incX, double const *Y, int const *incY);

    void dgemv_(char   const *TransA , int    const *M   , int    const *N  , double const *alpha ,
                double const *A      , int    const *lda , double const *X  , int    const *incX  ,
                double const *beta   , double       *Y   , int    const *incY);

    void dgemm_(char   const *TransA  , char   const *TransB , int    const *M    , int    const *N   ,
                int    const *K       , double const *alpha  , double const *A    , int    const *lda ,
                double const *B       , int    const *ldb    , double const *beta , double       *C   , int const *ldc);
}

namespace LinAlg
{


// LAPACK wrapper //////////////////////////////////////////////////////////////////////////////////////////////


/// internal product:  \f$ s \gets \{X\} \bullet \{Y\} \f$
inline double Dot(Vector<double> const & X, Vector<double> const & Y)
{
#ifndef DNDEBUG
    if (Y.Size()!=X.Size()) throw new Fatal("LinAlg::Axpy: The size (%d) of vector Y must be equal to the size (%d) of vector X",Y.Size(),X.Size());
#endif
    int N    = X.Size();
    int incX = 1;
    int incY = 1;
    return ddot_(&N,X.GetPtr(),&incX,Y.GetPtr(),&incY);
}

/// GEneral matrix * vector: \f$ \{Y\} \gets \alpha[A]\{X\} + \beta\{Y\}   \f$
inline void Gemv(double const a, Matrix<double> const & A, Vector<double> const & X, double const b, Vector<double> & Y)
{
    /*   {Y} <- a*[A]*{X} + b*{Y}
                    _                       _  / 0 \
         / 0 \     |  00 01 02 03 04 ... 0N  | | 1 |     / 0 \
         | 1 |     |  10 11 12 13 14 ... 1N  | | 2 |     | 1 |
         | 2 | = a*|  20 21 22 23 24 ... 2N  |*| 3 | + b*| 2 |
         |...|     |          ...            | | 4 |     |...|
         \ M /     |_ M0 M1 M2 M3 M4 ... MN _| |...|     \ M /
                                               \ N /
    
        This function preserves [A] matrix after multplication
    */
#ifndef DNDEBUG
    if (X.Size()!=A.Cols()) throw new Fatal("LinAlg::Gemv: The size (%d) of vector X must be equal to the number of columns (%d) of matrix A",X.Size(),A.Cols());
    if (Y.Size()!=A.Rows()) throw new Fatal("LinAlg::Gemv: The size (%d) of vector Y must be equal to the number of rows (%d) of matrix A",Y.Size(),A.Rows());
#endif
    int  M      = A.Rows();
    int  N      = A.Cols();
    char TransA = 'N';
    int  incX   = 1;
    int  incY   = 1;
    dgemv_(&TransA,        // CBLAS_TRANSPOSE A
           &M,             // M
           &N,             // N
           &a,             // Alpha
           A.GetPtr(),     // const double * A
           &M,             // LDA
           X.GetPtr(),     // const double * X
           &incX,          // incX
           &b,             // Beta
           Y.GetPtr(),     // double * Y
           &incY);         // incY
}

/// GEneral matrix x matrix: \f$ [C] \gets \alpha[A][B] + \beta[C]    \f$
inline void Gemm(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
/*  [C] <- a*[A]*[B] + b*[C] 
     _         _       _                       _  / 00 . 0N \      _         _
    |  00 . 0N  |     |  00 01 02 03 04 ... 0K  | | 10 . 1N |     |  00 . 0N  |
    |  10 . 1N  |     |  10 11 12 13 14 ... 1K  | | 20 . 2N |     |  10 . 1N  |
    |  20 . 2N  | = a*|  20 21 22 23 24 ... 2K  |*| 30 . 3N | + b*|  20 . 2N  |
    |    ...    |     |          ...            | | 40 . 4N |     |    ...    |
    |_ M0 . MN _|     |_ M0 M1 M2 M3 M4 ... MK _| |   ...   |     |_ M0 . MN _|
                                                  \ K0 . KN /
 */
#ifndef DNDEBUG
    if (B.Rows()!=A.Cols()) throw new Fatal("LinAlg::Gemm: The number of rows (%d) of matrix B must be equal to the number of columns (%d) of matrix A",B.Rows(),A.Cols());
    if (C.Rows()!=A.Rows()) throw new Fatal("LinAlg::Gemm: The number of rows (%d) of matrix C must be equal to the number of rows (%d) of matrix A",C.Rows(),A.Rows());
    if (C.Cols()!=B.Cols()) throw new Fatal("LinAlg::Gemm: The number of columns (%d) of matrix C must be equal to the number of columns (%d) of matrix B",C.Cols(),B.Cols());
#endif
    int  M      = A.Rows();
    int  K      = A.Cols();
    int  N      = B.Cols();
    char TransA = 'N';
    char TransB = 'N';
    dgemm_(&TransA,       // CBLAS_TRANSPOSE A
           &TransB,       // CBLAS_TRANSPOSE B
           &M,            // M
           &N,            // N
           &K,            // K
           &a,            // Alpha
           A.GetPtr(),    // const double * A
           &M,            // LDA
           B.GetPtr(),    // const double * B
           &K,            // LDB
           &b,            // Beta
           C.GetPtr(),    // double * C
           &M);           // LDC
}

/// GEneral matrix^T x matrix: \f$ [C] \gets \alpha[A]^T[B] + \beta[C]    \f$
inline void Gemtm(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
    if (B.Rows()!=A.Rows()) throw new Fatal("LinAlg::Gemtm: The number of rows (%d) of matrix B must be equal to the number of rows (%d) of matrix A",B.Rows(),A.Rows());
    if (C.Rows()!=A.Cols()) throw new Fatal("LinAlg::Gemtm: The number of rows (%d) of matrix C must be equal to the number of columns (%d) of matrix A",C.Rows(),A.Cols());
    if (C.Cols()!=B.Cols()) throw new Fatal("LinAlg::Gemtm: The number of columns (%d) of matrix C must be equal to the number of columns (%d) of matrix B",C.Cols(),B.Cols());
#endif
    int  K      = A.Rows();
    int  M      = A.Cols();
    int  N      = B.Cols();
    char TransA = 'T';
    char TransB = 'N';
    dgemm_(&TransA,       // CBLAS_TRANSPOSE A
           &TransB,       // CBLAS_TRANSPOSE B
           &M,            // M
           &N,            // N
           &K,            // K
           &a,            // Alpha
           A.GetPtr(),    // const double * A
           &K,            // LDA
           B.GetPtr(),    // const double * B
           &K,            // LDB
           &b,            // Beta
           C.GetPtr(),    // double * C
           &M);           // LDC
}

/// GEneral matrix x matrix^T: \f$ [C] \gets \alpha[A][B]^T+ \beta[C]    \f$
inline void Gemmt(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
    if (B.Cols()!=A.Cols()) throw new Fatal("LinAlg::Gemmt: The number of columns (%d) of matrix B must be equal to the number of columns (%d) of matrix A",B.Cols(),A.Cols());
    if (C.Rows()!=A.Rows()) throw new Fatal("LinAlg::Gemmt: The number of rows (%d) of matrix C must be equal to the number of rows (%d) of matrix A",C.Rows(),A.Rows());
    if (C.Cols()!=B.Rows()) throw new Fatal("LinAlg::Gemmt: The number of columns (%d) of matrix C must be equal to the number of rows (%d) of matrix B",C.Cols(),B.Rows());
#endif
    int  M      = A.Rows();
    int  K      = A.Cols();
    int  N      = B.Rows();
    char TransA = 'N';
    char TransB = 'T';
    dgemm_(&TransA,       // CBLAS_TRANSPOSE A
           &TransB,       // CBLAS_TRANSPOSE B
           &M,            // M
           &N,            // N
           &K,            // K
           &a,            // Alpha
           A.GetPtr(),    // const double * A
           &M,            // LDA
           B.GetPtr(),    // const double * B
           &N,            // LDB
           &b,            // Beta
           C.GetPtr(),    // double * C
           &M);           // LDC
}

/// GEneral matrix^T x matrix^T: \f$ [C] \gets \alpha[A]^T[B]^T + \beta[C]    \f$
inline void Gemtmt(double const a, Matrix<double> const & A, Matrix<double> const & B, double const b, Matrix<double> & C)
{
#ifndef DNDEBUG
    if (B.Cols()!=A.Rows()) throw new Fatal("LinAlg::Gemmt: The number of columns (%d) of matrix B must be equal to the number of rows (%d) of matrix A",B.Cols(),A.Rows());
    if (C.Rows()!=A.Cols()) throw new Fatal("LinAlg::Gemmt: The number of rows (%d) of matrix C must be equal to the number of columns (%d) of matrix A",C.Rows(),A.Cols());
    if (C.Cols()!=B.Rows()) throw new Fatal("LinAlg::Gemmt: The number of columns (%d) of matrix C must be equal to the number of rows (%d) of matrix B",C.Cols(),B.Rows());
#endif
    int  K      = A.Rows();
    int  M      = A.Cols();
    int  N      = B.Rows();
    char TransA = 'T';
    char TransB = 'T';
    dgemm_(&TransA,       // CBLAS_TRANSPOSE A
           &TransB,       // CBLAS_TRANSPOSE B
           &M,            // M
           &N,            // N
           &K,            // K
           &a,            // Alpha
           A.GetPtr(),    // const double * A
           &K,            // LDA
           B.GetPtr(),    // const double * B
           &N,            // LDB
           &b,            // Beta
           C.GetPtr(),    // double * C
           &M);           // LDC
}


// Implementation of basic operations //////////////////////////////////////////////////////////////////////////


// returns:  Y <- a*X
template <typename type>
inline
void _scale(size_t       n,
            type         a,
            type const * ptrX,
            type       * ptrY)
{
    size_t ll  = n % 4;
    for (size_t i = 0; i<ll; ++i) ptrY[i] = a*ptrX[i];
    for (size_t i = ll; i<n; i += 4)
    {
        ptrY[i]   = a*ptrX[i  ];
        ptrY[i+1] = a*ptrX[i+1];
        ptrY[i+2] = a*ptrX[i+2];
        ptrY[i+3] = a*ptrX[i+3];
    }
}

// returns:  Y <- a*X + b*W
template <typename type>
inline
void _scale(size_t       n,
            type         a,
            type const * ptrX,
            type         b,
            type const * ptrW,
            type       * ptrY)
{
    size_t ll  = n % 4;
    for (size_t i = 0; i<ll; ++i) ptrY[i] = a*ptrX[i] + b*ptrW[i];
    for (size_t i = ll; i<n; i += 4)
    {
        ptrY[i]   = a*ptrX[i  ] + b*ptrW[i  ];
        ptrY[i+1] = a*ptrX[i+1] + b*ptrW[i+1];
        ptrY[i+2] = a*ptrX[i+2] + b*ptrW[i+2];
        ptrY[i+3] = a*ptrX[i+3] + b*ptrW[i+3];
    }
}


// Scale function for Vector and Matrix ////////////////////////////////////////////////////////////////////////


// returns:  Y <- a*X
template <typename type>
inline
void scale(type                a,
          Vector<type> const & X,
          Vector<type>       & Y)
{
    size_t n = X.Size();
    Y.Resize(n);
    type const * ptrX = X.GetPtr();
    type       * ptrY = Y.GetPtr();
    _scale(n, a, ptrX, ptrY);
}
    
// returns:  Y <- a*X + b*W
template <typename type>
inline
void scale(type                a,
          Vector<type> const & X,
          type                 b,
          Vector<type> const & W,
          Vector<type>       & Y)
{
    int n = X.Size();
    if (W.Size()!=n) throw new Fatal("Internal Error: laexpr.h::scale: W.Size()==%d must be equal to %d",W.Size(),n);
    Y.Resize(n);
    type const * ptrX = X.GetPtr();
    type const * ptrW = W.GetPtr();
    type       * ptrY = Y.GetPtr();
    _scale(n, a, ptrX, b, ptrW, ptrY);
}

// returns:  Y <- a*X
template <typename type>
inline
void scale(type                a,
          Matrix<type> const & X,
          Matrix<type>       & Y)
{
    size_t n = X.Rows()*X.Cols();
    Y.Resize(X.Rows(), X.Cols());
    type const * ptrX = X.GetPtr();
    type       * ptrY = Y.GetPtr();
    _scale(n, a, ptrX, ptrY);
}
    
// returns:  Y <- a*X + b*W
template <typename type>
inline
void scale(type                a,
          Matrix<type> const & X,
          type                 b,
          Matrix<type> const & W,
          Matrix<type>       & Y)
{
    if (!(X.Rows()==W.Rows() && X.Cols()==W.Cols())) throw new Fatal("Internal Error: laexpr.h::scale: X.Rows==%d must be equal to W.Rows==%d and X.Cols==%d must be equal to W.Cols==%d",X.Rows(),W.Rows(),X.Cols(),W.Cols());
    size_t n = X.Rows()*X.Cols();
    Y.Resize(X.Rows(), X.Cols());
    type const * ptrX = X.GetPtr();
    type const * ptrW = W.GetPtr();
    type       * ptrY = Y.GetPtr();
    _scale(n, a, ptrX, b, ptrW, ptrY);
}

// returns:  B <- a*trn(A);
template <typename type>
inline
void scalet(type                 a,
            Vector<type> const & A,
            Matrix<type>       & B)
{
    size_t n = A.Size();
    B.Resize(1, n);
    type const * ptrA = A.GetPtr();
    type       * ptrB = B.GetPtr();
    _scale(n, a, ptrA, ptrB);
}

// returns:  B <- a*trn(A);
template <typename type>
inline
void scalet(type                 a,
            Matrix<type> const & A,
            Matrix<type>       & B) // n x m (col major)
{
    size_t m = A.Rows();
    size_t n = A.Cols();
    type const * ptrA;
    ptrA = A.GetPtr();

    B.Resize(n, m);
    type * ptrB = B.GetPtr();
    
    type rowA[n];
    for (size_t i=0; i<m; ++i)     //in rows of A
    {
        for (size_t k = 0; k<n; k++) //in cols of A
            rowA[k] = ptrA[i + k*m];
        type * colB = &ptrB[i*n];
        _scale(n, a, rowA, colB);
    }
}


// Function operations /////////////////////////////////////////////////////////////////////////////////////////


// Vector function operators

// Vector + Vector
template <typename type>
inline
Vector<type> & _add(Vector<type> const & A, Vector<type> const & B, Vector<type> & R)
{
    scale(static_cast<type>(1), A, static_cast<type>(1), B, R);
    return R;
}

// Vector - Vector
template <typename type>
inline
Vector<type> & _minus(Vector<type> const & A, Vector<type> const & B, Vector<type> & R)
{
    scale(static_cast<type>(1), A, static_cast<type>(-1), B, R);
    return R;
}

// scalar * Vector
template <class type, class t_num1>
inline
Vector<type> & _prod(t_num1 const & a, Vector<type> const & A, Vector<type> & R)
{
    //scale(static_cast<type>(a), A, R);
    R.Resize(A.Size());
    _scale(A.Size(), static_cast<type>(a), A.GetPtr(), R.GetPtr());
    return R;
}

// vector transpose 
template <typename type>
inline
Matrix<type> & _trn(Vector<type> const & A, Matrix<type> & R)
{
    scalet(static_cast<type>(1), A, R);
    return R;
}

// Matrix + Matrix
template <typename type>
inline
Matrix<type> & _add(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    scale(static_cast<type>(1), A, static_cast<type>(1), B, R);
    return R;
}

// Matrix - Matrix
template <typename type>
inline
Matrix<type> & _minus(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    scale(static_cast<type>(1), A, static_cast<type>(-1), B, R);
    return R;
}

// Matrix * Vector
template <class type>
inline
Vector<type> & _prod(Matrix<type> const & A, Vector<type> const & B, Vector<type> & R)
{
    R.Resize(A.Rows());
    if (A.Cols()!=B.Size()) throw new Fatal("Internal Error: laexpr.h::_prod: A.Cols==%d must be equal to B.Size==%d",A.Cols(),B.Size());
    Gemv(static_cast<type>(1), A, B, static_cast<type>(0), R);
    return R;
}

// Vector * Matrix
template <class type>
inline
Matrix<type> & _prod(Vector<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    // this operation is appropriate only when B is a row matrix
    size_t m = A.Size(); 
    size_t n = B.Cols();
    if (B.Rows()!=1) throw new Fatal("Internal Error: laexpr.h::_prod: B.Rows==%d must be equal to 1",B.Rows());
    R.Resize(m, n);
    type const * ptrA = A.GetPtr();
    type const * ptrB = B.GetPtr();
    type       * ptrR = R.GetPtr();

    for (size_t i=0; i<m; i++)
        for (size_t j=0; j<n; j++)
            ptrR[i+j*m] = ptrA[i]*ptrB[j];

    return R;
}

// Matrix * Matrix
template <class type>
inline
Matrix<type> & _prod(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    R.Resize(A.Rows(), B.Cols());
    Gemm(static_cast<type>(1), A, B, static_cast<type>(0), R);
    return R;
}

// scalar * Matrix
template <class type, class t_scalar>
inline
Matrix<type> & _prod(t_scalar const & A, Matrix<type> const & B, Matrix<type> & R)
{
    scale(static_cast<type const &>(A), B, R);
    return R;
}

// trn(Matrix)
template <typename type>
inline
Matrix<type> & _trn(Matrix<type> const & A, Matrix<type> & R)
{
    scalet(static_cast<type>(1), A, R);
    return R;
}

// trn(Matrix) * Matrix
template <class type>
inline
Matrix<type> & _prodt(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    R.Resize(A.Cols(), B.Cols());
    Gemtm(static_cast<type>(1), A, B, static_cast<type>(0), R);
    return R;
}

// Matrix * trn(Matrix)
template <class type>
inline
Matrix<type> & _prod_t(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    R.Resize(A.Rows(), B.Rows());
    Gemmt(static_cast<type>(1), A, B, static_cast<type>(0), R);
    return R;
}

// Vector * trn(Vector)
template <class type>
inline
Matrix<type> & _prod_t(Vector<type> const & A, Vector<type> const & B, Matrix<type> & R)
{
    size_t m = A.Size();
    size_t n = B.Size();
    R.Resize(m, n);
    type const * ptrA = A.GetPtr();
    type const * ptrB = B.GetPtr();
    type       * ptrR = R.GetPtr();
    for (size_t i=0; i<m; ++i)    
        for (size_t j=0; j<n; ++j)
            ptrR[i+j*m] = ptrA[i]*ptrB[j];
    return R;
}

// trn(Vector) * Vector
template <class type>
inline
type & _prodt(Vector<type> const & A, Vector<type> const & B, type & R)
{
    if (A.Size()!=B.Size()) throw new Fatal("Internal Error: laexpr.h::_prodt: A.Size==%d must be equal to B.Size==%d",A.Size(),B.Size());
    R = Dot(A, B);
    return R;
}

// trn(Matrix) * trn(Matrix)
template <class type>
inline
Matrix<type> & _prodtt(Matrix<type> const & A, Matrix<type> const & B, Matrix<type> & R)
{
    R.Resize(A.Cols(), B.Rows());
    Gemtmt(static_cast<type>(1), A, B, static_cast<type>(0), R);
    return R;
}


// Expression classes //////////////////////////////////////////////////////////////////////////////////////////


// Base expression class
template <class t_exp, class t_res>
class expression
{
public:
    operator t_res ()             const { return static_cast<t_exp const &>(*this).operator t_res(); }
    // Functions for external usage
    void Apply   (t_res & result) const { static_cast<t_exp const &>(*this).Apply   (result); }
    void Apply_pe(t_res & result) const { static_cast<t_exp const &>(*this).Apply_pe(result); }  
    void Apply_me(t_res & result) const { static_cast<t_exp const &>(*this).Apply_me(result); }  
protected:
    expression(){ }; // private constructor to avoid accidental constructions
};

// Output for expressions
template <class t_exp, class t_res>
std::ostream & operator << (std::ostream & os, expression<t_exp, t_res> const & expr)
{
    os << expr.operator t_res();
    return os;
}

// Expression class for unary operators
template <class t_exp1, class t_op, class t_res>
class exp_un:public expression<exp_un<t_exp1, t_op, t_res>, t_res>
{
public:
    typedef t_exp1  T_exp1;
    typedef t_op    T_op;
    typedef t_res   T_res;
    explicit exp_un(t_exp1 const & A):Arg(A){ } // explicit to avoid constructors ambiguity by default conversions
    operator t_res  () const { t_res result; return t_op::template Apply(Arg, result); }     
    void Apply(t_res & result) const 
    { 
        void const * ptr_arg    = static_cast<void const*>(&Arg);
        void const * ptr_result = static_cast<void const*>(&result);
        if (ptr_arg != ptr_result) t_op::template Apply(Arg, result); 
        else result = operator t_res();
    }    
    void Apply_pe(t_res & result) const { result = result + *this; }     
    void Apply_me(t_res & result) const { result = result - *this; }     
    t_exp1   const & Arg;
private:
    exp_un() { };  // private constructor to avoid accidental constructions
}; 

// Expression class for binary operators
template <class t_exp1, class t_exp2, class t_op, class t_res>
class exp_bin:public expression<exp_bin<t_exp1, t_exp2, t_op, t_res>, t_res>
{
public:
    typedef t_exp1 T_exp1;
    typedef t_exp2 T_exp2;
    typedef t_op   T_op;
    typedef t_res  T_res;
    exp_bin(t_exp1 const & A, t_exp2 const & B):Arg1(A),Arg2(B){}
    operator t_res () const { t_res result; return t_op::template Apply(Arg1, Arg2, result); }   
    void Apply(t_res & result) const 
    {
        void const * ptr_arg1   = static_cast<void const*>(&Arg1);
        void const * ptr_arg2   = static_cast<void const*>(&Arg2);
        void const * ptr_result = static_cast<void const*>(&result);
        if (ptr_arg1 != ptr_result && ptr_arg2 != ptr_result)
            t_op::template Apply(Arg1, Arg2, result); 
        else 
            result = operator t_res ();
    }    
    void Apply_pe(t_res & result) const { result = result + *this; }     
    void Apply_me(t_res & result) const { result = result - *this; }     
    t_exp1   const & Arg1;
    t_exp2   const & Arg2;
private:
    exp_bin() { };
};

// Expression class for ternary operators
template <class t_exp1, class t_exp2, class t_exp3, class t_op, class t_res>
class exp_ter:public expression<exp_ter<t_exp1, t_exp2, t_exp3, t_op, t_res>, t_res>
{
public:
    typedef t_exp1 T_exp1;
    typedef t_exp2 T_exp2;
    typedef t_exp3 T_exp3;
    typedef t_op   T_op;
    typedef t_res  T_res;
    exp_ter(t_exp1 const & A, t_exp2 const & B):Arg1(A),Arg2(B) { }
    operator t_res () const { t_res result; return t_op::template Apply(Arg1, Arg2, Arg3, result); }     
    void Apply(t_res & result) const 
    {
        void const * ptr_arg1   = static_cast<void const*>(&Arg1);
        void const * ptr_arg2   = static_cast<void const*>(&Arg2);
        void const * ptr_arg3   = static_cast<void const*>(&Arg3);
        void const * ptr_result = static_cast<void const*>(&result);
        if (ptr_arg1 != ptr_result && ptr_arg2 != ptr_result && ptr_arg3 != ptr_result)
            t_op::template Apply(Arg1, Arg2, Arg3, result); 
        else 
            result = operator t_res ();
    }    
    void Apply_pe(t_res & result) const { result = result + *this; }     
    void Apply_me(t_res & result) const { result = result - *this; }     
    t_exp1   const & Arg1;
    t_exp2   const & Arg2;
    t_exp3   const & Arg3;
private:
    exp_ter() { };
};


// Identification of expression result type ////////////////////////////////////////////////////////////////////


// for simple objects
template<class t_exp, class t_exp_again = t_exp>
struct res_type
{
    typedef t_exp T_res;
};
 
// for unary expressions
template<class t_exp>
struct res_type<t_exp, exp_un< typename t_exp::T_exp1, typename t_exp::T_op, typename t_exp::T_res> >
{
    typedef typename t_exp::T_res T_res;
};

// for binary expressions
template<class t_exp>
struct res_type<t_exp, exp_bin< typename t_exp::T_exp1, typename t_exp::T_exp2, typename t_exp::T_op, typename t_exp::T_res> >
{
    typedef typename t_exp::T_res T_res;
};


// Operator classes ////////////////////////////////////////////////////////////////////////////////////////////


// class operator for addition
class op_add
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    { 
        typedef typename res_type<t_exp1>::T_res t_res1;
        typedef typename res_type<t_exp2>::T_res t_res2;
        return _add( static_cast<t_res1 const &>(A), static_cast<t_res2 const &>(B), R); 
    }
};

// class operator for subtraction
class op_minus
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    {
        typedef typename res_type<t_exp1>::T_res t_res1;
        typedef typename res_type<t_exp2>::T_res t_res2;
        return _minus(static_cast<t_res1 const &>(A), static_cast<t_res2 const &>(B), R); 
    }
};

// class operator for unary minus
class op_minus_un
{
public:
    template <class t_exp1, class t_res>
    static t_res & Apply(t_exp1 const & A, t_res & R) 
    { 
        return _prod(-1.0, static_cast<t_res const &>(A), R); 
    }
};

// class operator for multiplication
class op_prod
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    { 
        typedef typename res_type<t_exp1>::T_res t_res1;
        typedef typename res_type<t_exp2>::T_res t_res2;
        return _prod( static_cast<t_res1 const &>(A), static_cast<t_res2 const &>(B), R); 
    }
};

class op_div_sc; //prototype for class operator that performs division by scalar

// class operator for multiplication by scalar 
class op_prod_sc 
{
public:
    template <class t_sc, class t_exp, class t_res>
    static t_res & Apply(t_sc const & A, t_exp const & B, t_res & R) 
    { 
        return _prod(A, static_cast<t_res const &>(B), R); 
    }

    // Overloading of apply functions to perform fast multiplication by multiple scalars:
    
    // scalar * exp_bin<op_prod_sc>
    template <class t_sc1, class t_sc2, class t_exp, class t_res>
    static t_res & Apply(t_sc1 const & A, exp_bin<t_sc2, t_exp, op_prod_sc, t_res > const & B, t_res & R) 
    { 
        return _prod(A*B.Arg1, static_cast<t_res const &>(B.Arg2), R); 
    }

    // scalar * exp_bin<op_div_sc>
    template <class t_sc1, class t_sc2, class t_exp, class t_res>
    static t_res & Apply(t_sc1 const & A, exp_bin<t_sc2, t_exp, op_div_sc, t_res > const & B, t_res & R) 
    { 
        return _prod(A/B.Arg1, static_cast<t_res const &>(B.Arg2), R); 
    }
};

// class operator for division by scalar 
class op_div_sc 
{
public:
    template <class t_sc, class t_exp, class t_res>
    static t_res & Apply(t_sc const & A, t_exp const & B, t_res & R) 
    { 
        return _prod(1.0/A, static_cast<t_res const &>(B), R); 
    }

    // Overloading of apply functions to perform fast multiplication by multiple scalars:
    
    // exp_bin<op_prod_sc> / scalar
    template <class t_sc1, class t_sc2, class t_exp, class t_res>
    static t_res & Apply(t_sc1 const & A, exp_bin<t_sc2, t_exp, op_prod_sc, t_res > const & B, t_res & R) 
    { 
        return _prod(B.Arg1/A, static_cast<t_res const &>(B.Arg2), R); 
    }

    // exp_bin<op_div_sc> / scalar
    template <class t_sc1, class t_sc2, class t_exp, class t_res>
    static t_res & Apply(t_sc1 const & A, exp_bin<t_sc2, t_exp, op_div_sc, t_res > const & B, t_res & R) 
    { 
        return _prod(1.0/(A*B.Arg1), static_cast<t_res const &>(B.Arg2), R); 
    }
};

// class operator for transposition
class op_trn
{
public:
    template <class t_exp1, class t_res>
    static t_res & Apply(t_exp1 const & A, t_res & R) 
    { 
        typedef typename res_type<t_exp1>::T_res t_res1;
        return _trn(static_cast<t_res1 const &>(A), R); 
    }
};

class op_oto // object transposed by object
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    { 
            //typedef typename res_type<t_exp1>::T_res t_res1;
            //typedef typename res_type<t_exp2>::T_res t_res2;
            _prodt(A, B, R);
            return R;
    }
};

class op_oot // object by object transposed
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    { 
            //typedef typename res_type<t_exp1>::T_res t_res1;
            //typedef typename res_type<t_exp2>::T_res t_res2;
            _prod_t(A, B, R);
            return R;
    }
};

class op_otot // object transposed by object transposed
{
public:
    template <class t_exp1, class t_exp2, class t_res>
    static t_res & Apply(t_exp1 const & A, t_exp2 const & B, t_res & R) 
    { 
            //typedef typename res_type<t_exp1>::T_res t_res1;
            //typedef typename res_type<t_exp2>::T_res t_res2;
            _prodtt(A, B, R);
            return R;
    }
};


// Expression & object operations //////////////////////////////////////////////////////////////////////////////


// macro for scalar operations
#define DECLARE_OP_WITH_SCALAR(MACRO) \
MACRO(short); \
MACRO(int); \
MACRO(long int); \
MACRO(size_t); \
MACRO(float); \
MACRO(double); \
MACRO(long double); 

// - expression (unary)
template <class t_exp, class t_res> 
inline  
exp_un< t_exp, op_minus_un, t_res >
operator-(expression<t_exp, t_res> const & A) 
{ 
  return exp_un< t_exp, op_minus_un, t_res >(static_cast<t_exp const &>(A)); 
} 

// + expression (unary)
template <class t_exp, class t_res> 
inline  
expression<t_exp, t_res>
operator+(expression<t_exp, t_res> const & A) 
{ 
  return A; 
} 

// expression + obj 
template <class t_obj, class t_exp1> 
inline  
exp_bin< t_exp1, t_obj, op_add, t_obj > 
operator+(expression<t_exp1, t_obj> const & A, t_obj const & B) 
{ 
  return exp_bin< t_exp1 , t_obj, op_add, t_obj >(static_cast<t_exp1 const &>(A),B); 
} 

// obj + expression
template <class t_obj, class t_exp1 > 
inline  
exp_bin< t_obj, t_exp1, op_add, t_obj >  
operator+(t_obj const & A, expression<t_exp1,t_obj> const & B) 
{ 
    return exp_bin<t_obj, t_exp1, op_add, t_obj >(A,static_cast<t_exp1 const &>(B)); 
} 

// exp + exp
template <class t_obj, class t_exp1, class t_exp2>
inline 
exp_bin< t_exp1, t_exp2, op_add, t_obj > 
operator+(expression<t_exp1, t_obj > const & A, expression<t_exp2, t_obj > const & B)
{
    return exp_bin<t_exp1 ,t_exp2, op_add, t_obj >(static_cast<t_exp1 const &>(A), static_cast<t_exp2 const &>(B));
}

// exp - exp
template <class t_obj, class t_exp1, class t_exp2>
inline 
exp_bin< t_exp1, t_exp2, op_minus, t_obj > 
operator-(expression<t_exp1, t_obj > const & A, expression<t_exp2, t_obj > const & B)
{
    return exp_bin<t_exp1 ,t_exp2, op_minus, t_obj >(static_cast<t_exp1 const &>(A), static_cast<t_exp2 const &>(B));
}

// expression - obj 
template <class t_obj, class t_exp1> 
inline  
exp_bin< t_exp1, t_obj, op_minus, t_obj > 
operator-(expression<t_exp1, t_obj> const & A, t_obj const & B) 
{ 
    return exp_bin< t_exp1 , t_obj, op_minus, t_obj >(static_cast<t_exp1 const &>(A),B); 
} 

// obj - expression
template <class t_obj, class t_exp1 > 
inline  
exp_bin< t_obj, t_exp1, op_minus, t_obj >  
operator-(t_obj const & A, expression<t_exp1,t_obj> const & B) 
{ 
    return exp_bin<t_obj, t_exp1, op_minus, t_obj >(A,static_cast<t_exp1 const &>(B)); 
} 

// scalar * expression && expression * scalar
#define EXPR_MULT_SCALAR(T)                                                           \
template <class t_obj, class t_exp1>                                                  \
inline exp_bin<T, t_exp1, op_prod_sc, t_obj >                                         \
operator*(T const & A, expression<t_exp1, t_obj > const & B)                          \
{                                                                                     \
    return exp_bin<T, t_exp1, op_prod_sc, t_obj >(A, static_cast<t_exp1 const &>(B)); \
}                                                                                     \
template <class t_obj, class t_exp1>                                                  \
inline exp_bin<T, t_exp1, op_prod_sc, t_obj >                                         \
operator*(expression<t_exp1, t_obj > const & A, T const & B)                          \
{                                                                                     \
    return exp_bin<T, t_exp1, op_prod_sc, t_obj >(B, static_cast<t_exp1 const &>(A)); \
} 

DECLARE_OP_WITH_SCALAR(EXPR_MULT_SCALAR);

// expression * scalar
#define EXPR_DIV_SCALAR(T)                                                            \
template <class t_obj, class t_exp1>                                                  \
inline exp_bin<T, t_exp1, op_div_sc, t_obj >                                          \
operator/(expression<t_exp1, t_obj > const & A, T const & B)                          \
{                                                                                     \
    return exp_bin<T, t_exp1, op_div_sc, t_obj >(B, static_cast<t_exp1 const &>(A));  \
} 

DECLARE_OP_WITH_SCALAR(EXPR_DIV_SCALAR);


// Vector operations ///////////////////////////////////////////////////////////////////////////////////////////


// - Vector
template <class type> 
inline  
exp_un< Vector<type>, op_minus_un, Vector<type> >
operator-(Vector<type> const & A) 
{ 
    return exp_un< Vector<type>, op_minus_un, Vector<type> > (A);
} 

// + Vector
template <class type> 
inline  
Vector<type>
operator+(Vector<type> const & A) 
{ 
    return A;
} 

// Vector + Vector
template <typename type>
inline 
exp_bin<Vector<type>, Vector<type>, op_add, Vector<type> >
operator + (Vector<type> const & A, Vector<type> const & B)
{
    return exp_bin<Vector<type>, Vector<type>, op_add, Vector<type> >(A,B);
}

// Vector - Vector
template <typename type>
inline 
exp_bin<Vector<type>, Vector<type>, op_minus, Vector<type> >
operator - (Vector<type> const & A, Vector<type> const & B)
{
    return exp_bin<Vector<type>, Vector<type>, op_minus, Vector<type> >(A,B);
}

// trn(Vector)
template <typename type>
inline
exp_un<Vector<type>, op_trn, Matrix<type> >
trn (Vector<type> const & A)
{
    return exp_un<Vector<type>, op_trn, Matrix<type> >(A);
}

// trn(exp)
template <class type, class t_exp1>
inline 
exp_un< t_exp1, op_trn, Matrix<type> > 
trn (expression<t_exp1, Vector<type> > const & A)
{
    return exp_un<t_exp1, op_trn, Matrix<type> >(static_cast<t_exp1 const &>(A));
}

//  scalar * Vector && Vector * scalar
#define VECTOR_MULT_SCALAR(T)                                           \
template <class type>                                                   \
inline                                                                  \
exp_bin< T, Vector<type>, op_prod_sc, Vector<type> >                    \
operator * ( T const & A, Vector<type> const & B)                       \
{                                                                       \
    return exp_bin< T, Vector<type>, op_prod_sc, Vector<type> >(A,B);   \
}                                                                       \
template <class type>                                                   \
inline                                                                  \
exp_bin<T, Vector<type>, op_prod_sc, Vector<type> >                     \
operator * (Vector<type> const & A, T const & B)                        \
{                                                                       \
    return exp_bin<T, Vector<type>, op_prod_sc, Vector<type> >(B,A);    \
}

DECLARE_OP_WITH_SCALAR(VECTOR_MULT_SCALAR);

//  Vector / scalar
#define VECTOR_DIV_SCALAR(T)                                            \
template <class type>                                                   \
inline                                                                  \
exp_bin<T, Vector<type>, op_div_sc, Vector<type> >                      \
operator / (Vector<type> const & A, T const & B)                        \
{                                                                       \
    return exp_bin<T, Vector<type>, op_div_sc, Vector<type> >(B,A);     \
}

DECLARE_OP_WITH_SCALAR(VECTOR_DIV_SCALAR);


// Matrix operations ///////////////////////////////////////////////////////////////////////////////////////////


// - Matrix
template <class type> 
inline  
exp_un< Matrix<type>, op_minus_un, Matrix<type> >
operator-(Matrix<type> const & A) 
{ 
    return exp_un< Matrix<type>, op_minus_un, Matrix<type> > (A);
} 

// + Matrix
template <class type> 
inline  
Matrix<type>
operator+(Matrix<type> const & A) 
{ 
    return A;
} 

// Matrix + Matrix
template <typename type>
inline 
exp_bin<Matrix<type>, Matrix<type>, op_add, Matrix<type> >
operator + (Matrix<type> const & A, Matrix<type> const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_add, Matrix<type> >(A,B);
}

// Matrix - Matrix
template <typename type>
inline 
exp_bin<Matrix<type>, Matrix<type>, op_minus, Matrix<type> >
operator - (Matrix<type> const & A, Matrix<type> const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_minus, Matrix<type> >(A,B);
}

// Matrix * Matrix
template <typename type>
inline 
exp_bin<Matrix<type>, Matrix<type>, op_prod, Matrix<type> >
operator * (Matrix<type> const & A, Matrix<type> const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_prod, Matrix<type> >(A,B);
}

// trn(Matrix)
template <typename type>
inline
exp_un<Matrix<type>, op_trn, Matrix<type> >
trn (Matrix<type> const & A)
{
    return exp_un<Matrix<type>, op_trn, Matrix<type> >(A);
}

// trn(exp)
template <class type, class t_exp1>
inline 
exp_un< t_exp1, op_trn, Matrix<type> > 
trn (expression<t_exp1, Matrix<type> > const & A)
{
    return exp_un<t_exp1, op_trn, Matrix<type> >(static_cast<t_exp1 const &>(A));
}

// trans(Matrix)
template <typename type>
inline
exp_un<Matrix<type>, op_trn, Matrix<type> >
trans (Matrix<type> const & A)
{
    return exp_un<Matrix<type>, op_trn, Matrix<type> >(A);
}

// trans(exp)
template <class type, class t_exp1>
inline 
exp_un< t_exp1, op_trn, Matrix<type> > 
trans (expression<t_exp1, Matrix<type> > const & A)
{
    return exp_un<t_exp1, op_trn, Matrix<type> >(static_cast<t_exp1 const &>(A));
}

//  scalar * Matrix && Matrix * scalar
#define MATRIX_MULT_SCALAR(T)                                           \
template <class type>                                                   \
inline exp_bin< T ,Matrix <type>, op_prod_sc,Matrix <type> >            \
operator * (T const & A, Matrix <type> const & B)                       \
{                                                                       \
    return exp_bin< T,Matrix <type>, op_prod_sc,Matrix <type> >(A,B);   \
}                                                                       \
template <class type>                                                   \
inline exp_bin<T, Matrix <type>, op_prod_sc,Matrix <type> >             \
operator * (Matrix <type> const & A, T const & B)                       \
{                                                                       \
    return exp_bin<T, Matrix <type>, op_prod_sc,Matrix <type> >(B,A);   \
}

DECLARE_OP_WITH_SCALAR(MATRIX_MULT_SCALAR);

//  Matrix / scalar
#define MATRIX_DIV_SCALAR(T)                                            \
template <class type>                                                   \
inline exp_bin<T, Matrix <type>, op_div_sc,Matrix <type> >              \
operator / (Matrix <type> const & A, T const & B)                       \
{                                                                       \
    return exp_bin<T, Matrix <type>, op_div_sc,Matrix <type> >(B,A);    \
}

DECLARE_OP_WITH_SCALAR(MATRIX_DIV_SCALAR);


// Matrix & Vector operations //////////////////////////////////////////////////////////////////////////////////


// Vector * Matrix
template <class type>
inline 
exp_bin<Vector<type>, Matrix<type>, op_prod, Matrix<type> >
operator * (Vector<type> const & A, Matrix<type> const & B)
{
    return exp_bin<Vector<type>, Matrix<type>, op_prod, Matrix<type> >(A,B);
}

// Matrix * Vector
template <typename type>
inline 
exp_bin<Matrix<type>, Vector<type>, op_prod, Vector<type> >
operator * (Matrix<type> const & A, Vector<type> const & B)
{
    return exp_bin<Matrix<type>, Vector<type>, op_prod, Vector<type> >(A,B);
}

// expression * Vector
template <class type, class t_exp1>
inline 
exp_bin< t_exp1, Vector<type>, op_prod, Vector<type> >
operator*(expression<t_exp1, Matrix<type> > const & A, Vector<type> const & B)
{
    return exp_bin< t_exp1 , Vector<type>, op_prod, Vector<type> >(static_cast<t_exp1 const &>(A),B);
}

// expression * Matrix
template <class type, class t_exp1>
inline 
exp_bin< t_exp1, Matrix<type>, op_prod, Matrix<type> >
operator*(expression<t_exp1, Matrix<type> > const & A, Matrix<type> const & B)
{
    return exp_bin< t_exp1, Matrix<type>, op_prod, Matrix<type> >(static_cast<t_exp1 const &>(A),B);
}

// Matrix * expression
template <class type, class t_exp1>
inline 
exp_bin<Matrix<type>, t_exp1, op_prod, Matrix<type> >
operator*(Matrix<type> const & A, expression<t_exp1, Matrix<type> > const & B)
{
    return exp_bin< Matrix<type>, t_exp1, op_prod, Matrix<type> >(A,static_cast<t_exp1 const &>(B));
}

// expression * expression (1)
template <class type, class t_exp1, class t_exp2>
inline 
exp_bin< t_exp1, t_exp2, op_prod, Matrix<type> >
operator*(expression<t_exp1, Matrix<type> > const & A, expression<t_exp2, Matrix<type> > const & B)
{
    return exp_bin< t_exp1, t_exp2, op_prod, Matrix<type> >(static_cast<t_exp1 const &>(A),static_cast<t_exp2 const &>(B));
}

// expression * expression (2)
template <class type, class t_exp1, class t_exp2>
inline 
exp_bin< t_exp1, t_exp2, op_prod, Vector<type> >
operator*(expression<t_exp1, Matrix<type> > const & A, expression<t_exp2, Vector<type> > const & B)
{
    return exp_bin< t_exp1, t_exp2, op_prod, Vector<type> >(static_cast<t_exp1 const &>(A),static_cast<t_exp2 const &>(B));
}

// expression * expression (3)
template <class type, class t_exp1, class t_exp2>
inline 
exp_bin< t_exp1, t_exp2, op_prod, Matrix<type> >
operator*(expression<t_exp1, Vector<type> > const & A, expression<t_exp2, Matrix<type> > const & B)
{
    return exp_bin< t_exp1, t_exp2, op_prod, Matrix<type> >(static_cast<t_exp1 const &>(A),static_cast<t_exp2 const &>(B));
}


// Matrix, Vector combinations with transposition //////////////////////////////////////////////////////////////


// trn(Vector) * Vector
template <class type>
inline 
exp_bin<Vector<type>, Vector<type>, op_oto, type >
operator * (exp_un< Vector<type>, op_trn, Matrix<type> > const & A, Vector<type> const & B)
{
    return exp_bin<Vector<type>, Vector<type>, op_oto, type >(A.Arg, B);
}

// Vector * trn(Vector)
template <class type>
inline 
exp_bin<Vector<type>,Vector<type>, op_oot, Matrix<type> >
operator * (Vector<type> const & A, exp_un< Vector<type>, op_trn, Matrix<type> > const & B)
{
    return exp_bin<Vector<type>,Vector<type>, op_oot, Matrix<type> >(A, B.Arg);
}

// trn(Matrix) * Matrix 
template <class type>
inline
exp_bin<Matrix<type>, Matrix<type>, op_oto, Matrix<type> >
operator * (exp_un<Matrix<type>, op_trn, Matrix<type> > const & A, Matrix<type> const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_oto, Matrix<type> >(A.Arg, B);
}

// Matrix * trn(Matrix)
template <class type>
inline
exp_bin<Matrix<type>, Matrix<type>, op_oot, Matrix<type> >
operator * (Matrix<type> const & A, exp_un<Matrix<type>, op_trn, Matrix<type> > const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_oot, Matrix<type> >(A, B.Arg);
}

// trn(Matrix) * trn(Matrix)
template <class type>
inline
exp_bin<Matrix<type>, Matrix<type>, op_otot, Matrix<type> >
operator * (exp_un<Matrix<type>, op_trn, Matrix<type> > const & A, exp_un<Matrix<type>, op_trn, Matrix<type> > const & B)
{
    return exp_bin<Matrix<type>, Matrix<type>, op_otot, Matrix<type> >(A.Arg, B.Arg);
}

}; // namespace LinAlg

#endif // MECHSYS_LINALG_LAEXPR_H
