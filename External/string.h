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

#ifndef MECHSYS_STRING_H
#define MECHSYS_STRING_H

// STL
#include <string>
#include <boost/algorithm/string.hpp> // for to_upper
#include <cstdarg> // for va_list, va_start, va_end
#include <cstdio>  // for vsnprintf
#include <stdio.h> // for printf

// color codes for the terminal
#define TERM_RST           "[0m"
#define TERM_CLR_BLACK     "[30m"
#define TERM_CLR_BLACK_H   "[1;30m"
#define TERM_CLR_WHITE     "[37m"
#define TERM_CLR_WHITE_H   "[1;37m"
#define TERM_CLR_RED       "[31m"
#define TERM_CLR_RED_H     "[1;31m"
#define TERM_CLR_GREEN     "[32m"
#define TERM_CLR_GREEN_H   "[1;32m"
#define TERM_CLR_YELLOW    "[33m"
#define TERM_CLR_YELLOW_H  "[1;33m"
#define TERM_CLR_BLUE      "[34m"
#define TERM_CLR_BLUE_H    "[1;34m"
#define TERM_CLR_MAGENTA   "[35m"
#define TERM_CLR_MAGENTA_H "[1;35m"
#define TERM_CLR_CYAN      "[36m"
#define TERM_CLR_CYAN_H    "[1;36m"
#define TERM_WHITE_BLUE    "[1;37;44m"
#define TERM_YELLOW_BLUE   "[1;33;44m"
#define TERM_YELLOW_BLACK  "[1;33;40m"
#define TERM_BLACK_WHITE   "[1;30;47m"

#ifdef TERM_WHITEBG
  #define TERM_CLR1  TERM_CLR_BLACK_H
  #define TERM_CLR2  TERM_CLR_BLACK
  #define TERM_CLR3  TERM_CLR_BLUE
  #define TERM_CLR4  TERM_CLR_CYAN
  #define TERM_CLR5  TERM_CLR_MAGENTA
  #define TERM_RED   TERM_CLR_RED_H
  #define TERM_GREEN TERM_CLR_GREEN
#else
  #ifdef TERM_NOCOLORS
    #define TERM_CLR1  ""
    #define TERM_CLR2  ""
    #define TERM_CLR3  ""
    #define TERM_CLR4  ""
    #define TERM_CLR5  ""
    #define TERM_RED   ""
    #define TERM_GREEN ""
  #else
    #define TERM_CLR1  TERM_CLR_YELLOW_H
    #define TERM_CLR2  TERM_CLR_WHITE_H
    #define TERM_CLR3  TERM_CLR_BLUE_H
    #define TERM_CLR4  TERM_CLR_CYAN_H
    #define TERM_CLR5  TERM_CLR_MAGENTA_H
    #define TERM_RED   TERM_CLR_RED_H
    #define TERM_GREEN TERM_CLR_GREEN_H
  #endif
#endif

class String : public std::string
{
public:
    // Constructors
    String () {}
    String (std::string const & Other) : std::string(Other) {}
    String (char        const * Other) : std::string(Other) {}

    // Methods
    int          Printf     (String const & Fmt, ...);                                                        ///< Print with format
    int          PrintfV    (String const & Fmt, va_list ArgList) { return _set_msg (Fmt.c_str(), ArgList); } ///< Print with format and ArgList
    char const * CStr       ()                              const { return this->c_str(); }                   ///< Get C-string
    void         TextFmt    (char const * NF);                                                                ///< Convert NF (ex: "%10g") to text Format (ex: "%10s")
    void         Split      (String & Left, String & Right, char const * Separator=" ");                      ///< Split string into left and right parts separated by Separator
    bool         HasWord    (String const & Word) { return (find(Word)!=npos); }                              ///< Check if string has a word Word
    void         GetFNKey   (String & FNKey);                                                                 ///< Return string without ending ".something"
    void         ToUpper    ()       { boost::to_upper ((*this)); }                                           ///< Convert string to upper case
    void         ToLower    ()       { boost::to_lower ((*this)); }                                           ///< Convert string to lower case
    String       ToUpperCpy () const { String tmp((*this)); tmp.ToUpper(); return tmp; }                      ///< (copy) Convert string to lower case
    String       ToLowerCpy () const { String tmp((*this)); tmp.ToLower(); return tmp; }                      ///< (copy) Convert string to lower case

    // For compatibility with wxWidgets
    String       & ToStdString()       { return (*this); }
    String const & ToStdString() const { return (*this); }

#ifdef USE_BOOST_PYTHON
    BPy::str PyStr() const { return BPy::str(CStr()); }
#endif

private:
    int _set_msg (char const * Fmt, va_list ArgList);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline int String::Printf (String const & Fmt, ...)
{
    int len;
    va_list       arg_list;
    va_start     (arg_list, Fmt);
    len=_set_msg (Fmt.c_str(), arg_list);
    va_end       (arg_list);
    return len;
}

inline void String::TextFmt (char const * NF)
{
    // number format for text
    this->clear();
    this->append(NF);
    size_t pos;
    pos=this->find("g"); while (pos!=String::npos) { this->replace(pos,1,"s"); pos=this->find("g",pos+1); }
    pos=this->find("f"); while (pos!=String::npos) { this->replace(pos,1,"s"); pos=this->find("f",pos+1); }
    pos=this->find("e"); while (pos!=String::npos) { this->replace(pos,1,"s"); pos=this->find("e",pos+1); }
}

inline void String::Split (String & Left, String & Right, char const * Separator)
{
    size_t pos = find(Separator);
    Left  = substr (0,pos);
    if (pos==npos) Right = "";
    else           Right = substr (pos+String(Separator).size());
}

inline void String::GetFNKey (String & FNKey)
{
    size_t pos = rfind(".");
    FNKey = substr (0,pos);
}

inline int String::_set_msg (char const * Fmt, va_list ArgList)
{
    const int size = 4048; // TODO: remove this limitation by using a loop and reallocating space
    char      buffer[size];
    int       len = std::vsnprintf(buffer, size, Fmt, ArgList);
    this->clear();
    if (len<0) this->append("String::_set_msg: INTERNAL ERROR: std::vsnprintf FAILED");
    else
    {
        buffer[len]='\0';
        if (len>size) this->append("String::_set_msg: INTERNAL ERROR: std::vsnprintf MESSAGE TRUNCATED: ");
        this->append(buffer);
    }
    return len;
}

#endif // MECHSYS_STRING_H
