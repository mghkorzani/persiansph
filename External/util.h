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

#ifndef MECHSYS_SORT_H
#define MECHSYS_SORT_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON
#include <fstream>
#include <sstream> // for istringstream

// MechSys
#include "string.h"
#include "array.h"

namespace Util
{

// Constants
const double ZERO   = sqrt(DBL_EPSILON); ///< Machine epsilon (smaller positive)
const double SQ2    = sqrt(2.0);         ///< \f$ \sqrt{2} \f$
const double SQ3    = sqrt(3.0);         ///< \f$ \sqrt{3} \f$
const double SQ6    = sqrt(6.0);         ///< \f$ \sqrt{6} \f$
const double SQ2BY3 = sqrt(2.0/3.0);     ///< \f$ \sqrt{2/3} \f$
const double PI     = 4.0*atan(1.0);     ///< \f$ \pi \f$

inline bool IsNan(double Val)
{
    return (std::isnan(Val) || ((Val==Val)==false)); // NaN is the only value, for which the expression Val==Val is always false
}

/** Signum function. */
inline double Signum (double x, double Tol=DBL_EPSILON) { return (fabs(x)>Tol ? (x>=0.0 ? +1.0 : -1.0) : 0.0); }

/*
inline double Sgn        (double Val)              { return (Val>=0.0 ? +1.0 : -1.0);                                } ///< Sgn function where Sgn(0)=+1
inline double Sign       (double a, double b)      { return (b>=0.0 ? fabs(a) : -fabs(a));                           } ///< Composite Sgn function. Returns |a| or -|a| according to the sign of b
inline double Acos       (double Val)              { return (Val>=1.0 ?  0.0 : (Val<=-1.0 ? Pi() : acos(Val)) );     } ///< Safe acos function
inline bool   Str2Bool   (String const & Str)      { if (Str==String("true") || Str==String("TRUE") || Str==String("True")) return true; else return false; } ///< Converts "TRUE", "True", or "true" to bool
inline bool   IsNanOrInf (double Val)              { int r=std::fpclassify(Val); if (r==FP_NAN) return true; if (r==FP_INFINITE) return true; return false; } ///< Check whether a number is NaN of Inf
*/

template <typename Val_T> inline Val_T Min (Val_T const & a, Val_T const & b) { return (a<b ? a : b); } ///< Minimum between a and b
template <typename Val_T> inline Val_T Max (Val_T const & a, Val_T const & b) { return (a>b ? a : b); } ///< Maximum between a and b
template <typename Val_T> inline Val_T Max (Val_T const & a, Val_T const & b, Val_T const & c)
{
    if (a>=b && a>=c) return a;
    if (b>=a && b>=c) return b;
    return c;
} ///< Maximum between a and b and c

/** Swap two values. */
template <typename Val_T>
inline void Swap (Val_T & a, Val_T & b)
{
    Val_T tmp = a;
    a = b;
    b = tmp;
}

/** Sort two values on an ascending order. */
template <typename Val_T>
inline void Sort (Val_T & a, Val_T & b)
{
    if (b<a) Util::Swap (a,b);
}

/** Sort three values on an ascending order. */
template <typename Val_T>
inline void Sort (Val_T & a, Val_T & b, Val_T & c)
{
    if (b<a) Util::Swap (a,b);
    if (c<b) Util::Swap (b,c);
    if (b<a) Util::Swap (a,b);
}

/** Sort four values on an ascending order. */
template <typename Val_T>
inline void Sort (Val_T & a, Val_T & b, Val_T & c, Val_T & d)
{
    if (b<a) Util::Swap (b,a);
    if (c<b) Util::Swap (c,b);
    if (d<c) Util::Swap (d,c);
    if (b<a) Util::Swap (b,a);
    if (c<b) Util::Swap (c,b);
    if (b<a) Util::Swap (b,a);
}

/** Find best square for given rows and columns. */
inline void FindBestSquare (int Size, int & nRow, int & nCol)
{
    nRow = -1;  // not found
    nCol = -1;  // not found
    for (int x=1; x<=Size; ++x)
    {
        if ((x*x)>=Size)
        {
            if ((x*x)==Size)
            {
                nRow = x;
                nCol = x;
                return;
            }
            else
            {
                for (int y=x; y>=1; --y)
                {
                    if ((x*y)==Size)
                    {
                        nRow = x;
                        nCol = y;
                        return;
                    }
                }
            }
        }
    }
}

struct FmtErr
{
    FmtErr (double TheError, double TheTol, char const * TheFmt="%g") : Error(TheError), Tol(TheTol), Fmt(TheFmt) {}
    double Error;
    double Tol;
    String Fmt;
};

std::ostream & operator<< (std::ostream & os, FmtErr const & P)
{
    String str;
    str.Printf (P.Fmt.CStr(), P.Error);
    os << (P.Error>P.Tol ? "[1;31m" : "[1;32m") << str << "[0m";
    return os;
}

inline bool FileExists (String const & Filename)
{
    std::ifstream file(Filename.CStr(), std::ios::in);
    if (file.fail()) return false;
    else
    {
        file.close();
        return true;
    }
}

inline bool HasKey (String const & KeysSepBySpace, String const & Key)
{
    std::istringstream iss(KeysSepBySpace);
    String key;
    bool has_key = false;
    while (iss>>key && !has_key) has_key = (key==Key);
    return has_key;
}

inline void Keys2Array (String const & KeysSepBySpace, Array<String> & Keys)
{
    std::istringstream iss(KeysSepBySpace);
    String key;
    while (iss>>key) Keys.Push (key);
}


#ifdef USE_BOOST_PYTHON

inline BPy::tuple PyFindBestSquare (int Size)
{
    int row, col;
    FindBestSquare (Size, row, col);
    return BPy::make_tuple (row, col);
}

inline void GetPyMethod (char const * ClassName, char const * MethodName, BPy::object & Method, char const * Filename="__main__")
{
    try
    {
        BPy::object main_module((BPy::handle<>(BPy::borrowed(PyImport_AddModule("__main__")))));
        BPy::object main_namespace = main_module.attr("__dict__");
        if (strcmp(Filename,"__main__")!=0)
        {
            if (!FileExists(Filename)) throw new Fatal("Util::GetPyMethod: Could not file named %s",Filename);
            BPy::exec_file (Filename, main_namespace, main_namespace);
        }
        BPy::object py_class        = BPy::extract<BPy::object>(main_namespace[ClassName])();
        BPy::object class_namespace = py_class.attr("__dict__");
        Method = BPy::extract<BPy::object>(class_namespace[MethodName])();
    }
    catch (BPy::error_already_set const & Err)
    {
        printf("\n%sUtil::GetPyMethod: Could not get method=='%s' of class=='%s' in '__main__'\nPython error message: ",TERM_CLR_RED_H,MethodName,ClassName);
        PyErr_Print();
        printf("%s\n",TERM_RST);
        throw new Fatal("PyODESolver::Init: failed (see message above).");
    }
}

inline void GetPyMethod (char const * InstanceName, char const * ClassName, char const * MethodName, BPy::object & Instance, BPy::object & Method, char const * Filename="__main__")
{
    try
    {
        BPy::object main_module((BPy::handle<>(BPy::borrowed(PyImport_AddModule("__main__")))));
        BPy::object main_namespace = main_module.attr("__dict__");
        if (strcmp(Filename,"__main__")!=0)
        {
            if (!FileExists(Filename)) throw new Fatal("Util::GetPyMethod: Could not file named %s",Filename);
            BPy::exec_file (Filename, main_namespace, main_namespace);
        }
        Instance = BPy::extract<BPy::object>(main_namespace[InstanceName])();
        BPy::object py_class        = BPy::extract<BPy::object>(main_namespace[ClassName])();
        BPy::object class_namespace = py_class.attr("__dict__");
        Method = BPy::extract<BPy::object>(class_namespace[MethodName])();
    }
    catch (BPy::error_already_set const & Err)
    {
        printf("\n%sUtil::GetPyMethod: Could not get method=='%s' of object=='%s' (class=='%s') in '__main__'\nPython error message: ",TERM_CLR_RED_H,MethodName,InstanceName,ClassName);
        PyErr_Print();
        printf("%s\n",TERM_RST);
        throw new Fatal("PyODESolver::Init: failed (see message above).");
    }
}

#endif


}; // namespace Util

#endif // MECHSYS_SORT_H
