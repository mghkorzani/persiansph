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

#ifndef MECHSYS_LINALG_VECTOR_H
#define MECHSYS_LINALG_VECTOR_H

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

/** Dense vector.
 Examples:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tmv.cpp?view=markup">   tmv.cpp Test matrix-vector multiplication</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsv.cpp?view=markup">   tsv.cpp Test LAPACK solver</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsymmv.cpp?view=markup">tsymmv.cpp Test symmetric matrix-vector multiplication</a>
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tvv.cpp?view=markup">   tvv.cpp Test vector-vector operations</a>
*/
template<typename Value_T>
class Vector
{
public:
    // Constructors
    Vector () : _size(0), _values(NULL), data(_values) {}  ///< Default constructor
    Vector (int Size);                                     ///< Constructor setting the size
    Vector (Vector<Value_T> const & Other);                ///< Copy constructor

    /** Destructor. */
    ~Vector () { if (_values!=NULL) delete [] _values; }

    // Access methods
    int             Size   () const; ///< Return the size of this vector
    Value_T const * GetPtr () const; ///< Return a pointer to the values, that cannot be modified, of this vector
    Value_T       * GetPtr ();       ///< Return a pointer to the values of this vector

    // Methods
    void Resize     (int Size);       ///< Resize this vector AND fill with ZERO values
    void change_dim (int Size) { Resize(Size); }
    void SetValues  (Value_T Value);  ///< Set all values equal to a given Value

    // Operators
    void            operator=  (Vector<Value_T> const & R); ///< Assignment operator
    void            operator+= (Vector<Value_T> const & R); ///< Plus-assignment operator
    void            operator-= (Vector<Value_T> const & R); ///< Minus-assignment operator
    void            operator/= (Value_T const & Scalar);    ///< Division operator
    void            operator*= (Value_T const & Scalar);    ///< Multiplication operator
    Value_T &       operator() (int i);                     ///< Write i value of vector
    const Value_T & operator() (int i) const;               ///< Read i value of vector

    // Auxiliar structures
    // operator to assign values separated by commas \cond
    class CommaAssign
    {
    public:
        CommaAssign(Value_T * Values, int & Size, Value_T const & FirstValue):_values(Values),_size(Size),_index(0) 
        { 
            if (_values==NULL) throw new Fatal("Vector::CommaAssign::CommaAssign (_values==NULL). The vector must be resized before calling this method.");
            _values[0] = FirstValue;
            for (int i=1; i<_size; ++i)
                _values[i] = static_cast<Value_T>(0);
        }
        CommaAssign & operator, (Value_T const & Num) 
        {
            _index++;
            if (_index>=_size) throw new Fatal("Vector::CommaAssign::operator, (_index>=_size). There are too many values in the comma expression to be assigned to this vector.");
            _values[_index] = Num;
            return *this;
        }
    private:
        Value_T * _values;
        int     & _size;
        int       _index;
    };

    CommaAssign operator= (Value_T const & Num) 
    {
        return CommaAssign(_values, _size, Num);
    }

    // Methods required by the expressions evaluation
    template<typename t_exp>
    void operator= (const expression<t_exp, Vector<Value_T> > & Exp) { Exp.Apply(*this); } 

    template<typename t_exp>
    void operator+= (const expression<t_exp, Vector<Value_T> > & Exp) { Exp.Apply_pe(*this); } 

    template<typename t_exp>
    void operator-= (const expression<t_exp, Vector<Value_T> > & Exp) { Exp.Apply_me(*this); }
    // \endcond

private:
    // Data
    int       _size;   ///< The number of components of this vector
    Value_T * _values; ///< The values of the components of this vector

public:
    Value_T * data;

}; // class Vector


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructors

template<typename Value_T>
inline Vector<Value_T>::Vector(int Size)
    : _size(0), _values(NULL), data(_values)
{
    if (Size>0)
    {
        _size   = Size;
        _values = new Value_T [_size];
        for (int i=0; i<_size; ++i) _values[i] = static_cast<Value_T>(0);
        data = _values;
    }
}

template<typename Value_T>
inline Vector<Value_T>::Vector(Vector<Value_T> const & R)
    : _size(R.Size()), _values(NULL), data(_values)
{
    if (_size>0)
    {
        // allocate memory
        _values = new Value_T [_size];
        data    = _values;

        // copy
        int i = 1;
        int j = 1;
        dcopy_ (&_size, R.GetPtr(), &i, _values, &j);
    }
}

// Access methods

template<typename Value_T>
inline int Vector<Value_T>::Size() const
{
    return _size;
}

template<typename Value_T>
inline const Value_T * Vector<Value_T>::GetPtr() const
{
#ifndef DNDEBUG
    //if (_values==NULL) throw new Fatal("Vector::Size: (_values==NULL). The vector must be resized before getting the pointer to its values.");
#endif
    return _values;
}

template<typename Value_T>
inline Value_T * Vector<Value_T>::GetPtr()
{
#ifndef DNDEBUG
    //if (_values==NULL) throw new Fatal("Vector::Size: (_values==NULL). The vector must be resized before getting the pointer to its values.");
#endif
    return _values;
}

// Methods
    
template<typename Value_T>
inline void Vector<Value_T>::Resize(int Size)
{
    if (Size>0)
    {
        if (Size==_size && _values!=NULL) return;   
        if (_values!=NULL) delete [] _values;
        _size   = Size;
        _values = new Value_T [_size];
        data    = _values;
        for (int i=0; i<_size; ++i) _values[i] = static_cast<Value_T>(0);
    }
    else
    {
        if (_values!=NULL) delete [] _values;
        _size   = 0;
        _values = NULL;
        data    = _values;
    }
}

template<typename Value_T>
inline void Vector<Value_T>::SetValues(Value_T Value)
{
#ifndef DNDEBUG
    if (_values==NULL) throw new Fatal("Vector::SetValues: (_values==NULL). The vector must be resized before setting its values.");
#endif
    for (int i=0; i<_size; ++i)
        _values[i] = Value; 
}

// Operators

template<typename Value_T>
inline void Vector<Value_T>::operator= (Vector<Value_T> const & R)
{
#ifndef DNDEBUG
    if (&R==this) throw new Fatal("Vector::operator= The right-hand-size of this operation (LHS = RHS) must not be equal to the LHS.");
#endif

    // reallocate if they are different (LHS != RHS)
    if (_values==NULL || R.Size()!=_size)
    {
        _size = R.Size();
        if (_values!=NULL) delete [] _values;
        _values = new Value_T [_size];
        data    = _values;
    }

    // copy
    int i = 1;
    int j = 1;
    dcopy_ (&_size, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Vector<Value_T>::operator+= (Vector<Value_T> const & R)
{
#ifndef DNDEBUG
    if (_values==NULL  ) return;//throw new Fatal("Vector::operator+= (_values==NULL). The vector must be resized before calling this method.");
    if (R.Size()!=_size) throw new Fatal("Vector::operator+= (R.Size()!=_size). The number of components of the LHS (%d) must be equal to the number of components of the RHS (%d).",R.Size(),_size);
#endif
    Value_T a = 1.0;
    int     i = 1;
    int     j = 1;
	daxpy_ (&_size, &a, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Vector<Value_T>::operator-= (Vector<Value_T> const & R)
{
#ifndef DNDEBUG
    if (_values==NULL  ) return;//throw new Fatal("Vector::operator-= (_values==NULL). The vector must be resized before calling this method.");
    if (R.Size()!=_size) throw new Fatal("Vector::operator-= (R.Size()!=_size). The number of components of the LHS (%d) must be equal to the number of components of the RHS (%d).",R.Size(),_size);
#endif
    Value_T a = -1.0;
    int     i = 1;
    int     j = 1;
	daxpy_ (&_size, &a, R.GetPtr(), &i, _values, &j);
}

template<typename Value_T>
inline void Vector<Value_T>::operator/= (Value_T const & Scalar)
{
#ifndef DNDEBUG
    if (_values==NULL) return;//throw new Fatal("Vector::operator/= (_values==NULL). The vector must be resized before calling this method.");
#endif
    Value_T a = 1.0/Scalar;
    int     d = 1;
    dscal_ (&_size, &a, _values, &d);
}

template<typename Value_T>
inline void Vector<Value_T>::operator*= (Value_T const & Scalar)
{
#ifndef DNDEBUG
    if (_values==NULL) return;//throw new Fatal("Vector::operator*= (_values==NULL). The vector must be resized before calling this method.");
#endif
    Value_T a = Scalar;
    int     d = 1;
    dscal_ (&_size, &a, _values, &d);
}

template<typename Value_T>
inline Value_T & Vector<Value_T>::operator() (int i)
{
#ifndef DNDEBUG
    if (_values==NULL  ) throw new Fatal("Vector::operator() (_values==NULL). The vector must be resized before calling this method.");
    if (i<0 || i>=_size) throw new Fatal("Vector::operator() (i<0 || i>=_size). The index (i=%d) for a component must be greater than/or equal to zero and smaller than the number of components (size=%d).",i,_size);
#endif
    return _values[i];
}

template<typename Value_T>
inline const Value_T & Vector<Value_T>::operator() (int i) const
{
#ifndef DNDEBUG
    if (_values==NULL  ) throw new Fatal("Vector::operator() (_values==NULL). The vector must be resized before calling this method.");
    if (i<0 || i>=_size) throw new Fatal("Vector::operator() (i<0 || i>=_size). The index (i=%d) for a component must be greater than/or equal to zero and smaller than the number of components (size=%d).",i,_size);
#endif
    return _values[i];
}

/** Outputs a vector. */
template<typename Value_T>
std::ostream & operator<< (std::ostream & os, const LinAlg::Vector<Value_T> & V)
{
    for (int i=0; i<V.Size(); ++i) os << V(i) << " ";
    os << std::endl;
    return os;
}


}; // namespace LinAlg

#endif // MECHSYS_LINALG_VECTOR_H
