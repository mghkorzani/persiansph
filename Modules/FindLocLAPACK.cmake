#####################################################################################
# MechSys - A C++ library to simulate Mechanical Systems                            #
# Copyright (C) 2010 Sergio Galindo                                                 #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

SET(LocLAPACK_LIBRARY_SEARCH_PATH
  $ENV{PKG}/lapack-3.5.0/lib
  $ENV{PKG}/lapack-3.5.0)

FIND_LIBRARY(LocLAPACK_BLAS     NAMES blas     PATHS ${LocLAPACK_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(LocLAPACK_LAPACK   NAMES lapack   PATHS ${LocLAPACK_LIBRARY_SEARCH_PATH})

SET(LocLAPACK_FOUND 1)
FOREACH(var LocLAPACK_BLAS LocLAPACK_LAPACK)
  IF(NOT ${var})
	SET(LocLAPACK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(LocLAPACK_FOUND)
  SET(LocLAPACK_LIBRARIES ${LocLAPACK_LAPACK} ${LocLAPACK_BLAS})
ENDIF(LocLAPACK_FOUND)
