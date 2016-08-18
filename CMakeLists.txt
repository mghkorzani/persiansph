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

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

PROJECT                (New_Test)

INCLUDE ($ENV{SPH}/Modules/FindPKG.cmake)

SET(EXES
#	Drag_LES
#	Poiseuille
#	Overtopping
#	Hydraulic_Heave
#	Current-flow_elastoplastic
	Jet_Scouring
)

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)


