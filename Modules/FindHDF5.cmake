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

SET(HDF5_INCLUDE_SEARCH_PATH
    $ENV{PKG}/hdf5-1.8.16/hl/src
    $ENV{PKG}/hdf5-1.8.16/src
    /usr/include)

SET(HDF5_LIBRARY_SEARCH_PATH
    $ENV{PKG}/hdf5-1.8.16/hl/src/.libs
    $ENV{PKG}/hdf5-1.8.16/src/.libs
    /usr/lib/x86_64-linux-gnu
    /usr/lib)

FIND_PATH(HDF5_H    hdf5.h    ${HDF5_INCLUDE_SEARCH_PATH})
FIND_PATH(HDF5_HL_H hdf5_hl.h ${HDF5_INCLUDE_SEARCH_PATH})

#FIND_LIBRARY(Z       NAMES z       PATHS ${HDF5_LIBRARY_SEARCH_PATH})
#FIND_LIBRARY(SZ      NAMES szip    PATHS ${HDF5_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(HDF5    NAMES hdf5    PATHS ${HDF5_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(HDF5_HL NAMES hdf5_hl PATHS ${HDF5_LIBRARY_SEARCH_PATH})

SET(HDF5_FOUND 1)
#FOREACH(var HDF5_H HDF5_HL_H Z SZ HDF5 HDF5_HL)
FOREACH(var HDF5_H HDF5_HL_H HDF5 HDF5_HL)
  IF(NOT ${var})
	SET(HDF5_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(HDF5_FOUND)
  SET(HDF5_INCLUDE_DIR  ${HDF5_H} ${HDF5_HL_H})
#  SET(HDF5_LIBRARIES    ${HDF5_HL} ${HDF5} ${SZ} ${Z})
  SET(HDF5_LIBRARIES    ${HDF5_HL} ${HDF5})
ENDIF(HDF5_FOUND)

