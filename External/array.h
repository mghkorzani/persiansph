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

#ifndef MECHSYS_ARRAY_H
#define MECHSYS_ARRAY_H

// STL
#include <algorithm> // for std::find, std::min_element, and std::max_element
#ifdef USE_STDVECTOR
  #include <vector> // for std::vector
  #define INIT_VALUES
#else
  #define INIT_VALUES SzFactor=1.2; _values=NULL;
#endif

// MechSys
#include "fatal.h"

template<typename Value_T>
class Array
#ifdef USE_STDVECTOR
    : public std::vector<Value_T>
#endif
{
public:
    // Alternative constructors
    Array (Value_T const & v0, bool JustOne);
    Array (Value_T const & v0, Value_T const & v1);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18, Value_T const & v19);

#ifdef USE_BOOST_PYTHON
    // Python constructor
    Array (BPy::list const & Dat);
#endif

#ifdef USE_STDVECTOR
    // Constructors
    Array () {}                           ///< Default constructor
    Array (size_t Size) { Resize(Size); } ///< Alternative constructor

    // Methods
    size_t          Size   () const { return  std::vector<Value_T>::size();  } ///< Returns the size
    Value_T       * GetPtr ()       { return &std::vector<Value_T>::front(); } ///< Returns a pointer to the values
    Value_T const * GetPtr () const { return &std::vector<Value_T>::front(); } ///< Returns a pointer to the values
    Value_T       & Last   ()       { return  std::vector<Value_T>::back();  } ///< Return the last element
    Value_T const & Last   () const { return  std::vector<Value_T>::back();  } ///< Return the last element

    // NOTE:
    //  std::vector operators [] do not properly check for wrong indices when accessing the values.
    //  This sometimes causes segmentation fault.
    //  To solve this we could wrap operator[] and check the range of indices, but this does not work,
    //  since it fails to return a reference when using std::vector<bool> 

#else
    // Constructors and destructor
     Array ()            { INIT_VALUES Resize(0); }     ///< Default constructor
     Array (size_t Size) { INIT_VALUES Resize(Size); }  ///< Alternative constructor
     Array (Array<Value_T> const & Other);              ///< Copy constructor (needed when using Array< Array<...> >)
    ~Array () { if (_values!=NULL) delete [] _values; } ///< Destructor

    // Methods
    size_t          Size   () const { return _size; }            ///< Returns the size
    Value_T       * GetPtr ()       { return _values; }          ///< Returns a pointer to the values
    Value_T const * GetPtr () const { return _values; }          ///< Returns a pointer to the values
    Value_T       & Last   ()       { return _values[_size-1]; } ///< Return the last element
    Value_T const & Last   () const { return _values[_size-1]; } ///< Return the last element
    size_t          size   () const { return Size(); }           ///< Alternative Size() method

    // Operators
    Value_T       & operator[] (size_t i);                 ///< Access operator (write)
    Value_T const & operator[] (size_t i) const;           ///< Access operator (read)
    void            operator=  (Array<Value_T> const & R); ///< Assignment operator (needed when using Array< Array<...> >)
#endif

    // Methods
    void            Resize    (size_t Size);                       ///< Resize the array
    void            Push      (Value_T const & Value);             ///< Add a new entry increasing the size if necessary
    void            XPush     (Value_T const & Value);             ///< Exclusive Push: push only if Value is not already in array (not fast since it calls the Has method)
    void            PushN     (Value_T const & Value, size_t Num); ///< Add a new entry increasing the size if necessary
    long            Find      (Value_T const & Value) const;       ///< Find a value: returns -1 if not found, otherwise, returns the index of the element found
    bool            Has       (Value_T const & Value) const;       ///< Has Value ~ Find(Value)>=0 ?
    Value_T const & TheMin    () const;                            ///< Find the minimum value
    Value_T const & TheMax    () const;                            ///< Find the maximum value
    Value_T         Mean      () const;                            ///< Calculate the mean value (Value_T must have addition operators)
    Value_T         Norm      () const;                            ///< Calculate the norm value (Value_T must have addition operators)
    void            SetValues (Value_T const & V);                 ///< Set all values to be equal to V
    void            Clear     () { Resize(0); }                    ///< Clear array
    void            DelItem   (size_t i);                          ///< Delete item i from array (not efficient)
    void            DelItems  (Array<int> const & Idxs);           ///< Delete items from array (not efficient)
    void            DelVal    (Value_T const & Value);             ///< Delete Value from array (not efficient)

    // Assign values separated by commas
    class CommaAssign
    {
    public:
        CommaAssign (Array<Value_T> * ptArray, Value_T const & FirstValue) : _ptarray(ptArray), _i(0)
        {
            if (_ptarray->Size()>_i) (*_ptarray)[_i] = FirstValue;
            else throw new Fatal("Array::CommaAssign: The array must be resized with a size greater than or equal to %d before calling this method.",_i+1);
        }
        CommaAssign & operator, (Value_T const & Value)
        {
            _i++;
            if (_ptarray->Size()>_i) (*_ptarray)[_i] = Value;
            else throw new Fatal("Array::CommaAssign: The array must be resized with a size greater than or equal to %d before calling this method.",_i+1);
            return *this;
        }
    private:
        Array<Value_T> * _ptarray;
        size_t           _i;
    };
    CommaAssign operator= (Value_T const & Value) { return CommaAssign(this, Value); }

#ifndef USE_STDVECTOR
    double    SzFactor; ///< Scaling factor
private:
    size_t    _size;    ///< Current number of components
    size_t    _space;   ///< Available space
    Value_T * _values;  ///< Space to hold all values
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifndef USE_STDVECTOR

template<typename Value_T>
inline Array<Value_T>::Array(Array<Value_T> const & Other)
{
    INIT_VALUES
    Resize(Other.Size());
    for (size_t i=0; i<_size; ++i)
        _values[i] = Other[i];
}

#endif

#ifdef USE_BOOST_PYTHON

template<typename Value_T>
inline Array<Value_T>::Array (BPy::list const & Dat)
{
    INIT_VALUES
    size_t size = BPy::len(Dat);
    Resize (size);
    for (size_t i=0; i<size; ++i)
        (*this)[i] = BPy::extract<Value_T>(Dat[i])();
}

#endif

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, bool JustOne)
{
    INIT_VALUES
    Resize(1);
    (*this)[0] = v0;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1)
{
    INIT_VALUES
    Resize(2);
    (*this)[0] = v0;
    (*this)[1] = v1;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2)
{
    INIT_VALUES
    Resize(3);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3)
{
    INIT_VALUES
    Resize(4);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4)
{
    INIT_VALUES
    Resize(5);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5)
{
    INIT_VALUES
    Resize(6);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
    (*this)[5] = v5;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6)
{
    INIT_VALUES
    Resize(7);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
    (*this)[5] = v5;
    (*this)[6] = v6;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7)
{
    INIT_VALUES
    Resize(8);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
    (*this)[5] = v5;
    (*this)[6] = v6;
    (*this)[7] = v7;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8)
{
    INIT_VALUES
    Resize(9);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
    (*this)[5] = v5;
    (*this)[6] = v6;
    (*this)[7] = v7;
    (*this)[8] = v8;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9)
{
    INIT_VALUES
    Resize(10);
    (*this)[0] = v0;
    (*this)[1] = v1;
    (*this)[2] = v2;
    (*this)[3] = v3;
    (*this)[4] = v4;
    (*this)[5] = v5;
    (*this)[6] = v6;
    (*this)[7] = v7;
    (*this)[8] = v8;
    (*this)[9] = v9;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10)
{
    INIT_VALUES
    Resize(11);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11)
{
    INIT_VALUES
    Resize(12);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12)
{
    INIT_VALUES
    Resize(13);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13)
{
    INIT_VALUES
    Resize(14);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14)
{
    INIT_VALUES
    Resize(15);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15)
{
    INIT_VALUES
    Resize(16);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
    (*this)[15] = v15;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16)
{
    INIT_VALUES
    Resize(17);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
    (*this)[15] = v15;
    (*this)[16] = v16;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17)
{
    INIT_VALUES
    Resize(18);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
    (*this)[15] = v15;
    (*this)[16] = v16;
    (*this)[17] = v17;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18)
{
    INIT_VALUES
    Resize(19);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
    (*this)[15] = v15;
    (*this)[16] = v16;
    (*this)[17] = v17;
    (*this)[18] = v18;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18, Value_T const & v19)
{
    INIT_VALUES
    Resize(20);
    (*this)[ 0] = v0;
    (*this)[ 1] = v1;
    (*this)[ 2] = v2;
    (*this)[ 3] = v3;
    (*this)[ 4] = v4;
    (*this)[ 5] = v5;
    (*this)[ 6] = v6;
    (*this)[ 7] = v7;
    (*this)[ 8] = v8;
    (*this)[ 9] = v9;
    (*this)[10] = v10;
    (*this)[11] = v11;
    (*this)[12] = v12;
    (*this)[13] = v13;
    (*this)[14] = v14;
    (*this)[15] = v15;
    (*this)[16] = v16;
    (*this)[17] = v17;
    (*this)[18] = v18;
    (*this)[19] = v19;
}


// Methods

template<typename Value_T>
inline void Array<Value_T>::Resize (size_t Size)
{
#ifdef USE_STDVECTOR
    std::vector<Value_T>::resize (Size);
#else
    // Check
    if (SzFactor<1.0) throw new Fatal("Array::Resize: SzFactor==%f must be greater than 1.0", SzFactor);

    // Clear previous memory
    if (_values!=NULL) delete [] _values;

    // Allocate new memory
    _size   = Size;
    _space  = static_cast<size_t>((_size+1)*SzFactor+1);
    _values = new Value_T [_space];
#endif
}

template<typename Value_T>
inline void Array<Value_T>::Push (Value_T const & Value)
{
#ifdef USE_STDVECTOR
    std::vector<Value_T>::push_back (Value);
#else
    if (_size==_space)
    {
        size_t oldsz = _size;
        Value_T * tmp = new Value_T [oldsz];
        for (size_t i=0; i<oldsz; ++i) tmp[i] = _values[i];
        Resize (oldsz+1);
        for (size_t i=0; i<oldsz; ++i) _values[i] = tmp[i];
        delete [] tmp;
    }
    else _size++;
    _values[_size-1] = Value;
#endif
}

template<typename Value_T>
inline void Array<Value_T>::XPush (Value_T const & Value)
{
    if (!Has(Value)) Push (Value);
}

template<typename Value_T>
inline void Array<Value_T>::PushN (Value_T const & Value, size_t Num)
{
#ifdef USE_STDVECTOR
    std::vector<Value_T>::reserve (Num);
    for (size_t i=0; i<Num; ++i) std::vector<Value_T>::push_back (Value);
#else
    if (_size+Num>=_space)
    {
        size_t oldsz = _size;
        Value_T * tmp = new Value_T [oldsz];
        for (size_t i=0; i<oldsz; ++i) tmp[i] = _values[i];
        Resize (oldsz+Num);
        for (size_t i=0; i<oldsz; ++i) _values[i] = tmp[i];
        delete [] tmp;
    }
    else _size += Num;
    for (size_t i=0; i<Num; ++i) _values[_size-Num+i] = Value;
#endif
}

template<typename Value_T>
inline long Array<Value_T>::Find (Value_T const & Value) const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::find (std::vector<Value_T>::begin(), std::vector<Value_T>::end(), Value);
    if (it==std::vector<Value_T>::end()) return -1;
    else return it-std::vector<Value_T>::begin();
#else
    Value_T * res = std::find(_values, _values+_size, Value);
    if (res==_values+_size) return -1;
    else return res-_values;
#endif
}

template<typename Value_T>
inline bool Array<Value_T>::Has (Value_T const & Value) const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::find (std::vector<Value_T>::begin(), std::vector<Value_T>::end(), Value);
    if (it==std::vector<Value_T>::end()) return false;
    else return true;
#else
    Value_T * res = std::find(_values, _values+_size, Value);
    if (res==_values+_size) return false;
    else return true;
#endif
}

template<typename Value_T>
inline Value_T const & Array<Value_T>::TheMin () const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::min_element (std::vector<Value_T>::begin(), std::vector<Value_T>::end());
    return (*it);
#else
    Value_T * res = std::min_element(_values, _values+_size);
    return (*res);
#endif
}

template<typename Value_T>
inline Value_T const & Array<Value_T>::TheMax () const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::max_element (std::vector<Value_T>::begin(), std::vector<Value_T>::end());
    return (*it);
#else
    Value_T * res = std::max_element(_values, _values+_size);
    return (*res);
#endif
}

/*
template<typename Value_T>
inline long Array<Value_T>::Min () const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::min_element (std::vector<Value_T>::begin(), std::vector<Value_T>::end());
    return it-std::vector<Value_T>::begin();
#else
    Value_T * res = std::min_element(_values, _values+_size);
    return res-_values;
#endif
}

template<typename Value_T>
inline long Array<Value_T>::Max () const
{
#ifdef USE_STDVECTOR
    typename std::vector<Value_T>::const_iterator it = std::max_element (std::vector<Value_T>::begin(), std::vector<Value_T>::end());
    return it-std::vector<Value_T>::begin();
#else
    Value_T * res = std::max_element(_values, _values+_size);
    return res-_values;
#endif
}
*/

template<typename Value_T>
inline Value_T Array<Value_T>::Mean () const
{
    Value_T sum = 0.0;
    for (size_t i=0; i<Size(); ++i) sum += (*this)[i];
    return sum/Size();
}

template<typename Value_T>
inline Value_T Array<Value_T>::Norm () const
{
    Value_T sum = 0.0;
    for (size_t i=0; i<Size(); ++i) sum += (*this)[i]*(*this)[i];
    return sqrt(sum);
}

template<typename Value_T>
inline void Array<Value_T>::SetValues (Value_T const & V)
{
    for (size_t i=0; i<Size(); ++i) (*this)[i] = V;
}

template<typename Value_T>
inline void Array<Value_T>::DelItem (size_t IdxToRemove)
{
#ifndef NDEBUG
    if (IdxToRemove<0 || IdxToRemove>=Size()) throw new Fatal("Array::DelItem: Subscript==%zd (size==%zd) is out of range.", IdxToRemove, Size());
#endif
    Array<Value_T> tmp((*this));
    this->Resize (tmp.Size()-1);
    size_t k = 0;
    for (size_t i=0; i<tmp.Size(); ++i)
    {
        if (i!=IdxToRemove)
        {
            (*this)[k] = tmp[i];
            k++;
        }
    }
}

template<typename Value_T>
inline void Array<Value_T>::DelItems (Array<int> const & Idxs)
{
    if (Idxs.Size()>Size()) throw new Fatal("Array::DelItems: Array 'Idxs' with indices to be deleted must have size smaller than or equal to %zd. Idxs.Size()==%zd is invalid.",Size(),Idxs.Size());
    Array<Value_T> tmp((*this));
    this->Resize (tmp.Size()-Idxs.Size());
    size_t k = 0;
    for (size_t i=0; i<tmp.Size(); ++i)
    {
        if (!Idxs.Has(i))
        {
            (*this)[k] = tmp[i];
            k++;
        }
    }
}

template<typename Value_T>
inline void Array<Value_T>::DelVal (Value_T const & V)
{
    size_t idx_to_remove = Find(V);
    if (idx_to_remove<0) throw new Fatal("Array::DelVal: Value cannot be found for deletion");
    DelItem (idx_to_remove);
}


// Operators

#ifndef USE_STDVECTOR

template<typename Value_T>
inline Value_T & Array<Value_T>::operator[] (size_t i)
{
#ifndef NDEBUG
    if (i<0 || i>=Size()) throw new Fatal("Array::operator[] (write) Subscript==%zd (size==%zd) is out of range.", i, Size());
#endif
    return _values[i];
}

template<typename Value_T>
inline Value_T const & Array<Value_T>::operator[] (size_t i) const
{
#ifndef NDEBUG
    if (i<0 || i>=Size()) throw new Fatal("Array::operator[] (read) Subscript==%zd (size==%zd) is out of range.", i, Size()); 
#endif
    return _values[i];
}

template<typename Value_T>
inline void Array<Value_T>::operator= (Array<Value_T> const & R)
{
#ifndef DNDEBUG
    if (&R==this) throw new Fatal("Array::operator= The right-hand-size of this operation (LHS = RHS) must not be equal to the LHS.");
#endif
    // Reallocate if they are different (LHS != RHS)
    if (_size!=R.Size()) Resize(R.Size());
    
    // Copy values
    for (size_t i=0; i<_size; ++i) _values[i] = R[i];
}

#endif


/** Outputs an array. */
template<typename Value_T>
std::ostream & operator<< (std::ostream & os, const Array<Value_T> & V)
{
    for (size_t i=0; i<V.Size(); ++i)
    {
        os << V[i];
        if (i!=V.Size()-1) os << ", ";
    }
    return os;
}

#undef INIT_VALUES

#endif // MECHSYS_ARRAY_H
