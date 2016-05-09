#####################################################################################
# PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids         # 
#             and soils) using Smoothed Particle Hydrodynamics method               #   
# Copyright (C) 2016 Maziar Gholami Korzani and Sergio Galindo-Torres               #
#                                                                                   #
# This file is part of PersianSPH                                                   #
#                                                                                   #
# This is free software; you can redistribute it and/or modify it under the         #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 3 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# PersianSPH; if not, see <http://www.gnu.org/licenses/>                            #
#####################################################################################

SET(LocGSL_INCLUDE_SEARCH_PATH
  $ENV{PKG}/gsl-2.1
  $ENV{PKG}/gsl-2.1/gsl
  $ENV{PKG}/gsl-2.1/eigen
  /usr/include
  /usr/local/include)

SET(LocGSL_LIBRARY_SEARCH_PATH
  $ENV{PKG}/gsl-2.1/.libs
  $ENV{PKG}/gsl-2.1/cblas/.libs
  /usr/lib
  /usr/local/lib)

FIND_PATH(GSL_MATH_H  gsl/gsl_math.h  ${LocGSL_INCLUDE_SEARCH_PATH})
FIND_PATH(GSL_EIGEN_H gsl/gsl_eigen.h ${LocGSL_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(GSL_GSL      NAMES gsl      PATHS ${LocGSL_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(GSL_CBLASGSL NAMES gslcblas PATHS ${LocGSL_LIBRARY_SEARCH_PATH})


SET(LocGSL_FOUND 1)
FOREACH(var GSL_MATH_H GSL_EIGEN_H GSL_GSL GSL_CBLASGSL)
  IF(NOT ${var})
	SET(LocGSL_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(GSL_FOUND)
  SET(LocGSL_INCLUDE_DIRS ${GSL_MATH_H} ${GSL_EIGEN_H})
  SET(LocGSL_LIBRARIES    ${GSL_GSL} ${GSL_CBLASGSL})
ENDIF(GSL_FOUND)
