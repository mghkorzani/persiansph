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

#ifndef SPH_SPECIAL_FUNCTIONS_H
#define SPH_SPECIAL_FUNCTIONS_H

#include "matvec.h"

namespace SPH {

	double Kernel(size_t Dim, size_t KT, double r, double h);

	double GradKernel(size_t Dim, size_t KT, double r, double h);

	double LaplaceKernel(size_t Dim, size_t KT, double r, double h);

	double SecDerivativeKernel(size_t Dim, size_t KT, double r, double h);

	double EOS(size_t EQ, double Cs0, double P00, double Density, double Density0);

	double SoundSpeed(size_t EQ, double Cs0, double Density, double Density0);

	double DensitySolid(size_t EQ, double Cs0, double P00, double Pressure, double Density0);

	void   Seepage(size_t ST, double k, double k2, double mu,  double rho, double& SF1, double& SF2);

	void   Rotation (Mat3_t Input, Mat3_t & Vectors, Mat3_t & VectorsT, Mat3_t & Values);

	Mat3_t abab (const Mat3_t & A, const Mat3_t & B);

}; // namespace SPH

#include "Functions.cpp"

#endif // SPH_SPECIAL_FUNCTIONS_H
