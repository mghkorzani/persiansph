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

#ifndef MECHSYS_FATAL_H
#define MECHSYS_FATAL_H

// Std lib
#include <iostream> // for cout
#include <cstdarg>  // for va_list, va_start, va_end
#include <exception>

// MechSys
#include "string.h"

// Global variable to tell 'catch' that parallel code is activated
bool MECHSYS_CATCH_PARALLEL = false;

// MPI
#ifdef HAS_MPI
  #include <mpi.h>
  #define MECHSYS_FINALIZE  if (MECHSYS_CATCH_PARALLEL) MPI::COMM_WORLD.Abort(666); return 1;
  #define MECHSYS_MPI_CATCH catch (MPI::Exception e) { printf("%sFatal: MPI Error # %d : %s%s\n",TERM_RED,e.Get_error_code(),e.Get_error_string(),TERM_RST); MECHSYS_FINALIZE }
  #define MECHSYS_MPI_INIT  MPI::Init(argc, argv); \
                            MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
#else
  #define MECHSYS_FINALIZE  return 1;
  #define MECHSYS_MPI_CATCH
  #define MECHSYS_MPI_INIT
#endif

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python/errors.hpp>
  #include <boost/python/str.hpp>
  namespace BPy = boost::python;
  #define MECHSYS_BPY_CATCH catch (BPy::error_already_set) { printf("%sFatal: ",TERM_RED); PyErr_Print(); printf("%s\n",TERM_RST); MECHSYS_FINALIZE }
#else
  #define MECHSYS_BPY_CATCH
#endif

// Catch structure
#define MECHSYS_CATCH catch (Fatal      * e)     { e->Cout();  delete e;                                                   MECHSYS_FINALIZE } \
                      catch (char const * m)     { printf("%sFatal: %s%s\n",TERM_RED,m       ,TERM_RST);                   MECHSYS_FINALIZE } \
                      catch (std::exception & e) { printf("%sFatal: %s%s\n",TERM_RED,e.what(),TERM_RST);                   MECHSYS_FINALIZE } \
                      MECHSYS_MPI_CATCH                                                                                                       \
                      MECHSYS_BPY_CATCH                                                                                                       \
                      catch (...)                { printf("%sFatal: Some exception (...) occurred%s\n",TERM_RED,TERM_RST); MECHSYS_FINALIZE }

class Fatal
{
public:
	// Constructor
	Fatal (String const & Fmt, ...);

	// Methods
	void   Cout () const { printf("%sFatal: %s%s\n", TERM_RED, _msg.CStr(), TERM_RST); }
	String Msg  () const { return _msg; }

private:
	String _msg;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Fatal::Fatal(String const & Fmt, ...)
{
	va_list       arg_list;
	va_start     (arg_list, Fmt);
	_msg.PrintfV (Fmt, arg_list);
	va_end       (arg_list);
}


#ifdef USE_BOOST_PYTHON

void PyExceptTranslator (Fatal * e)
{
	PyErr_SetString(PyExc_UserWarning, e->Msg().CStr());
	delete e;
}

#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FATAL_H
