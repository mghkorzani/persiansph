/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *
*             and soils) using Smoothed Particle Hydrodynamics method              *
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Boundary_Condition.h"

namespace SPH {

  inline Boundary::Boundary()
  {
    allv	= 0.0;
    inv		= 0.0;
    outv	= 0.0;
    Periodic[0]=Periodic[1]=Periodic[2] = false;
    inDensity = 0.0;
    outDensity = 0.0;
    allDensity = 0.0;
    InOutFlow = 0;
    InFlowLoc1 = 0.0;
    InFlowLoc2 = 0.0;
    InFlowLoc3 = 0.0;
    OutFlowLoc = 0.0;
    cellfac = 3.9;
    inoutcounter = 0;
    MassConservation = false;
  }

}; // namespace SPH
