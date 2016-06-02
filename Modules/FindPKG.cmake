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

OPTION(A_MAKE_VERBOSE       "Show additional messages during compilation/linking?" ON )
OPTION(A_MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                       ON )
OPTION(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"                         OFF)
OPTION(A_MAKE_OPTIMIZED     "Make optimized (-O3)"                                 ON )

SET (FLAGS   "${FLAGS}")
SET (LIBS     ${LIBS})
SET (LFLAGS  "${LFLAGS}")
SET (MISSING "")

IF(A_MAKE_VERBOSE)
	SET (CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(A_MAKE_VERBOSE)

IF(A_MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS (-Wall)
ENDIF(A_MAKE_ALL_WARNINGS)

IF(A_MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS (-g)
ENDIF(A_MAKE_DEBUG_SYMBOLS)

IF(A_MAKE_OPTIMIZED)
	ADD_DEFINITIONS (-O3)
ENDIF(A_MAKE_OPTIMIZED)

ADD_DEFINITIONS(-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.
#ADD_DEFINITIONS(-std=c++11)         # New C++ standard
ADD_DEFINITIONS(-fpermissive)       # New C++ standard

INCLUDE      ($ENV{SPH}/Modules/FindHDF5.cmake)
INCLUDE      (FindOpenMP)
INCLUDE      (FindLAPACK)

set(GSL_GLOB_PATH ${CMAKE_ROOT}/Modules)
FIND_PATH(GSL_GLOB  FindGSL.cmake ${GSL_GLOB_PATH})
IF(NOT ${GSL_Glob} MATCHES "NOTFOUND" )
	INCLUDE      (FindGSL)
ENDIF(NOT ${GSL_Glob} MATCHES "NOTFOUND")

INCLUDE_DIRECTORIES($ENV{SPH}/Source $ENV{SPH}/External $ENV{PKG}/blitz-0.9)

if(HDF5_FOUND)
        ADD_DEFINITIONS (-DH5_NO_DEPRECATED_SYMBOLS -DH5Gcreate_vers=2 -DH5Gopen_vers=2 -DUSE_HDF5)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
else(HDF5_FOUND)
        SET (MISSING "${MISSING} HDF5")
endif(HDF5_FOUND)


if(OPENMP_FOUND)
    ADD_DEFINITIONS (-DUSE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
        SET (MISSING "${MISSING} OpenMP")
endif(OPENMP_FOUND)


if(LAPACK_FOUND)
    SET (LIBS ${LIBS} ${LAPACK_LIBRARIES})
else(LAPACK_FOUND)
     INCLUDE ($ENV{SPH}/Modules/FindLocLAPACK.cmake)
    if(LocLAPACK_FOUND)
        SET (LIBS ${LIBS} ${LocLAPACK_LIBRARIES} "gfortran")
    else(LocLAPACK_FOUND)
        SET (MISSING "${MISSING} Lapack & Blas")
    endif(LocLAPACK_FOUND)
endif(LAPACK_FOUND)


if(GSL_FOUND)
	INCLUDE_DIRECTORIES (${GSL_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${GSL_LIBRARIES})
else(GSL_FOUND)
     INCLUDE ($ENV{SPH}/Modules/FindLocGSL.cmake)
    if(LocGSL_FOUND)
	INCLUDE_DIRECTORIES (${LocGSL_INCLUDE_DIRS})
        SET (LIBS ${LIBS} ${LocGSL_LIBRARIES})
    else(LocGSL_FOUND)
        SET (MISSING "${MISSING} Gsl & GslcBlas")
    endif(LocGSL_FOUND)
endif(GSL_FOUND)


if(MISSING)
    MESSAGE("!!!!!!!!!!!!!!!!!!!!!!!! Missing dependencies =${MISSING} !!!!!!!!!!!!!!!!!!!!!!!!!!!!")
endif(MISSING)

