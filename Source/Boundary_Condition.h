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

#ifndef SPH_BOUNDARY_CONDITION_H
#define SPH_BOUNDARY_CONDITION_H

namespace SPH {

  class Boundary
  {
  public:
    // Data
    double	inDensity;	///< Apply a certain density to inflow particles
    double	outDensity;	///< Apply a certain density to outflow particles
    double	allDensity;	///< Apply a certain density to outflow particles

    Vec3_t 	inv;		///< Apply a certain velocity to inflow particles
    Vec3_t 	outv;		///< Apply a certain velocity to outflow particles
    Vec3_t	allv;		///< Apply a certain velocity to all particle

    bool 	Periodic[3];	///< Considering periodic in all directions => 0=X, 1=Y, 2=Z

    int 	InOutFlow;	///< Considering inflow in all directions  by adding and deleting particles=> [0]=X, [1]=Y, [2]=Z and 0=none, 1=-
    double	InFlowLoc1;
    double	InFlowLoc2;
    double	InFlowLoc3;
    double	OutFlowLoc;
    double	cellfac;
    int		inoutcounter;
    bool	MassConservation;

    Array <int>	OutPart;
    Array <int>	InPart;

    Boundary();
  };
}; // namespace SPH

#include "Boundary_Condition.cpp"

#endif // SPH_BOUNDARY_CONDITION_H
