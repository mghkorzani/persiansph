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

#ifndef MECHSYS_LINALG_MATRIX_H
#define MECHSYS_LINALG_MATRIX_H

// STL
#include <iostream>
#include <cmath>

// MechSys
#include "fatal.h"

extern "C"
{
    // BLAS
    void dscal_(int const *N, double const *alpha, double *X, int const *incX);
	void dcopy_(int const *N, double const *X, int const *incX, double *Y, int const *incY);
	void daxpy_(int const *N, double const *alpha, double const *X, int const *incX, double *Y, int const *incY);
}

namespace LinAlg
{

// Prototype for expression classes
template <class t_exp, class t_res>
class expression; // should not be called (directly) by the user

/** Dense matrix.

   Note:
          The values are internally saved in Column-Major order,
          thus, the method GetPtr() may be passed to a FORTRAN-style
          routine with no problems.

 Examples:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tmv.cpp?view=markup">   tmv.cpp Test matrix-vector multiplication</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsv.cpp?view=markup">   tsv.cpp Test LAPACK solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsymmv.cpp?view=markup">tsymmv.cpp Test symmetric matrix-vector multiplication</a>
*/
template<typename Value_T>
class Matrix
{
public:
    // Constructors
    Matrix () : _rows(0), _cols(0), _values(NULL), data(_values) {} ///< Default constructor
    Matrix (int Rows, int Cols);                                    ///< Constructor setting Rows and Cols
    Matrix (Matrix<Value_T> const & Other);                         ///< Copy constructor

    /** Destructor. */
    ~Matrix () { if (_values!=NULL) delete [] _values; }

    // Access methods
    int             Rows     () const; ///< Returns the number of rows of this matrix
    int             Cols     () const; ///< Returns the number of columns of this matrix
    size_t          num_rows () const { return (size_t)Rows(); }
    size_t          num_cols () const { return (size_t)Cols(); }
    Value_T const * GetPtr   () const; ///< Returns a pointer to the values, that cannot be modified, of this matrix
    Value_T       * GetPtr   ();       ///< Returns a pointer to the values of this matrix

    // Methods
    void Resize     (int Rows, int Cols); ///< Resize this matrix AND fill with ZERO values
    void change_dim (int Rows, int Cols) { Resize(Rows,Cols); }
    void SetValues  (Value_T Value);      ///< Set all values equal to a given Value

    // Operators
    void            operator=  (Matrix<Value_T> const & R); ///< Assignment operator
    void            operator+= (Matrix<Value_T> const & R); ///< Plus-assignment operator
    void            operator-= (Matrix<Value_T> const & R); ///< Minus-assignment operator
    void            operator/= (Value_T const & Scalar);    ///< Division operator
    void            operator*= (Value_T const & Scalar);    ///< Multiplication operator
    Value_T &       operator() (int i, int j);              ///< Write ij value of matrix
    const Value_T & operator() (int i, int j) const;        ///< Read ij value of matrix
    
    // Auxiliar structures
    // operator to assign values separated by commas \cond
    class CommaAssign
    {
    public:
        CommaAssign(Value_T * Values, int & Rows, int & Cols, Value_T const & FirstValue):_values(Values),_rows(Rows),_cols(Cols),_index(0) 
        { 
            if (_values==NULL) throw new Fatal("Matrix::CommaAssign::CommaAssign (_values==NULL). The matrix must be resized prior to use this method.");
            _values[0] = FirstValue;
            int _components = _rows*_cols;
            for (int i=1; i<_components; ++i)
                _values[i] = static_cast<Value_T>(0);
        }
        CommaAssign & operator, (Value_T const & Num) 
        {
            _index++;
            if (_index>=_rows*_cols) throw new Fatal("Matrix::CommaAssign::operator, (_index>=_rows*_cols). There are too many values in the comma expression to be assigned to this matrix.");
            _values[_index/_cols + (_index % _cols) * _rows] = Num; // considering col major
            return *this;
        }
    private:
        Value_T * _values;
        int     & _rows;
        int     & _cols;
        int       _index;
    };

    CommaAssign operator= (Value_T const & Num) 
    {
        return CommaAssign(_values, _rows, _cols, Num);
    }

    // Methods required by the expressions evaluation
    template<typename t_exp>
    void operator= (const expression<t_exp, Matrix<Value_T> > & Exp) { Exp.Apply(*this); } 

    template<typename t_exp>
    void operator+= (const expression<t_exp, Matrix<Value_T> > & Exp) { Exp.Apply_pe(*this); } 

    template<typename t_exp>
    void operator-= (const expression<t_exp, Matrix<Value_T> > & Exp) { Exp.Apply_me(*this); }
    // \endcond

private:
    // Data
    int       _rows;   ///< Number of rows
    int       _cols;   ///< Number of columns
    Value_T * _values; ///< Values in column-major order (like FORTRAN)

public:
    Value_T * data;

}; // class Matrix


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructors

template<typename Value_T>
inline Matrix<Value_T>::Matrix(int Rows, int Cols)
    : _rows(0), _cols(0), _values(NULL), data(_values)
{
    if (Rows>0 && Cols>0)
    {
        _rows = Rows;
        _cols = Cols;
        _values = new Value_T [_rows*_cols];
        for (int i=0; i<_rows*_cols; ++i)
            _values[i] = static_cast<Value_T>(0);
        data = _values;
    }
}

template<typename Value_T>
inline Matrix<Value_T>::Matrix(Matrix<Value_T> const & R)
    : _rows(R.Rows()), _cols(R.Cols()), _values(NULL), data(_values)
{
    if (_rows>0 && _cols>0)
    {
        // allocate memory
        _values = new Value_T [_rows*_cols];
        data    = _values;

        // copy
        int n = _rows*_cols;
        int i = 1;
        int j = 1;
        dcopy_ (&n, R.GetPtr(), &i, _values, &j);
    }
}

// Access methods

template<typename Value_T>
inline int Matrix<Value_T>::Rows() const
{
    return _rows;
}

template<typename Value_T>
inline int Matrix<Value_T>::Cols() const
{
    return _cols;
}

template<typename Value_T>
inline Value_T * Matrix<Value_T>::GetPtr()
{
#ifndef DNDEBUG
    //if (_values==NULL) throw new Fatal("Matrix::GetPtr: (_values==NULL). The matrix must be resized prior to get a pointer to its values.");
#endif
    return _values;
}

template<typename Value_T>
inline const Value_T * Matrix<Value_T>::GetPtr() const
{
#ifndef DNDEBUG
    //if (_values==NULL) throw new Fatal("Matrix::GetPtr: (_values==NULL). The matrix must be resized prior to get a pointer to its values.");
#endif
    return _values;
}

// Methods

template<typename Value_T>
inline void Matrix<Value_T>::Resize(int Rows, int Cols)
{
    if (Rows>0 && Cols>0)
    {
        if (Rows==_rows && Cols==_cols && _values!=NULL) return;    
        if (_values!=NULL) delete [] _values;
        _rows   = Rows;
        _cols   = Cols;
        _values = new Value_T [_rows*_cols];
        data    = _values;
        for (int i=0; i<_rows*_cols; ++i) _values[i] = static_cast<Value_T>(0);
    }
    else
    {
        if (_values!=NULL) delete [] _values;
        _rows   = 0;
        _cols   = 0;
        _values = NULL;
        data    = _values;
    }
}

template<typename Value_T>
inline void Matrix<Value_T>::SetValues(Value_T Value)
{
#ifndef DNDEBUG
    if (_values==NULL) throw new Fatal("Matrix::SetValues: (_values==NULL). The matrix must be resized prior to set its values.");
#endif
    for (int i=0; i<_rows*_cols; ++i)
        _values[i] = Value; 
}

// Operators

template<typename Value_T>
inline void Matrix<Value_T>::operator= (Matrix<Value_T> const & R)
{
#ifndef DNDEBUG
    if (&R==this) throw new Fatal("Matrix::operator= The right-hand-size of this operation (LHS = RHS) must not be equal to the LHS.");
#endif

    // reallocate if they are different (LHS != RHS)
    if (_values==NULL || R.Rows()!=_rows || R.Cols()!=_cols)
    {
        _rows = R.Rows();
        _cols = R.Cols();
        if (_values!=NULL) delete [] _values;
        _values = new Value_T [_rows*_cols];
        data    = _values;
    }

    // copy
    int n = _rows*_cols;
    int i = 1;
    int j = 1;
    dcopy_ (&n, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Matrix<Value_T>::operator+= (Matrix<Value_T> const & R)
{
#ifndef DNDEBUG
    if (_values==NULL  ) return;//throw new Fatal("Matrix::operator+= (_values==NULL). The matrix must be resized prior to use this method.");
    if (R.Rows()!=_rows) throw new Fatal("Matrix::operator+= (R.Rows()!=_rows). The number of Rows of the LHS (%d) must be equal to the number of rows of the RHS (%d).",R.Rows(),_rows);
    if (R.Cols()!=_cols) throw new Fatal("Matrix::operator+= (R.Rows()!=_cols). The number of Cols of the LHS (%d) must be equal to the number of cols of the RHS (%d).",R.Cols(),_cols);
#endif
    int     n = _rows*_cols;
    Value_T a = 1.0;
    int     i = 1;
    int     j = 1;
	daxpy_ (&n, &a, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Matrix<Value_T>::operator-= (Matrix<Value_T> const & R)
{
#ifndef DNDEBUG
    if (_values==NULL  ) return;//throw new Fatal("Matrix::operator-= (_values==NULL). The matrix must be resized prior to use this method.");
    if (R.Rows()!=_rows) throw new Fatal("Matrix::operator-= (R.Rows()!=_rows). The number of Rows of the LHS (%d) must be equal to the number of rows of the RHS (%d).",R.Rows(),_rows);
    if (R.Cols()!=_cols) throw new Fatal("Matrix::operator-= (R.Rows()!=_cols). The number of Cols of the LHS (%d) must be equal to the number of cols of the RHS (%d).",R.Cols(),_cols);
#endif
    int     n = _rows*_cols;
    Value_T a = -1.0;
    int     i = 1;
    int     j = 1;
	daxpy_ (&n, &a, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Matrix<Value_T>::operator/= (Value_T const & Scalar)
{
#ifndef DNDEBUG
    if (_values==NULL) return;//throw new Fatal("Matrix::operator/= (_values==NULL). The matrix must be resized before calling this method.");
#endif
    int     n = _rows*_cols;
    Value_T a = 1.0/Scalar;
    int     d = 1;
    dscal_ (&n, &a, _values, &d);
}

template<typename Value_T>
inline void Matrix<Value_T>::operator*= (Value_T const & Scalar)
{
#ifndef DNDEBUG
    if (_values==NULL) return;//throw new Fatal("Matrix::operator*= (_values==NULL). The matrix must be resized before calling this method.");
#endif
    int     n = _rows*_cols;
    Value_T a = Scalar;
    int     d = 1;
    dscal_ (&n, &a, _values, &d);
}

template<typename Value_T>
inline Value_T & Matrix<Value_T>::operator() (int i, int j)
{
#ifndef DNDEBUG
    if (_values==NULL  ) throw new Fatal("Matrix::operator() (_values==NULL). The matrix must be resized prior to use this method.");
    if (i<0 || i>=_rows) throw new Fatal("Matrix::operator() (i<0 || i>=_rows). The index (i=%d) for a row must be greater than/or equal to zero and smaller than the number of rows (nrow=%d).",i,_rows);
    if (j<0 || j>=_cols) throw new Fatal("Matrix::operator() (j<0 || j>=_cols). The index (j=%d) for a column must be greater than/or equal to zero and smaller than the number of columns (ncol=%d).",j,_cols);
#endif
    return _values[j*_rows+i]; // Col major
}

template<typename Value_T>
inline const Value_T & Matrix<Value_T>::operator() (int i, int j) const
{
#ifndef DNDEBUG
    if (_values==NULL  ) throw new Fatal("Matrix::(const)operator() (_values==NULL). The matrix must be resized prior to use this method.");
    if (i<0 || i>=_rows) throw new Fatal("Matrix::(const)operator() (i<0 || i>=_rows). The index (i=%d) for a row must be greater than/or equal to zero and smaller than the number of rows (nrow=%d).",i,_rows);
    if (j<0 || j>=_cols) throw new Fatal("Matrix::(const)operator() (j<0 || j>=_cols). The index (j=%d) for a column must be greater than/or equal to zero and smaller than the number of columns (ncol=%d).",j,_cols);
#endif
    return _values[j*_rows+i]; // Col major
}

/** Outputs a dense matrix. */
template<typename Value_T>
std::ostream & operator<< (std::ostream & os, const Matrix<Value_T> & M)
{
    for (int i=0; i<M.Rows(); ++i)
    {
        for (int j=0; j<M.Cols(); ++j) os << M(i,j) << " ";
        os << std::endl;
    }
    return os;
}


}; // namespace LinAlg

#endif // MECHSYS_LINALG_MATRIX_H
