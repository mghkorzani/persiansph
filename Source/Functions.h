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

	double Kernel								(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	double GradKernel						(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	double LaplaceKernel				(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	double SecDerivativeKernel	(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	double EOS									(size_t const & EQ, double const & Cs0, double const & P00, double const & Density, double const & Density0);

	double SoundSpeed						(size_t const & EQ, double const & Cs0, double const & Density, double const & Density0);

	double DensitySolid					(size_t const & EQ, double const & Cs0, double const & P00, double const & Pressure, double const & Density0);

	void   Seepage							(size_t const & ST, double const & k, double const & k2, double const & mu,  double const & rho, double & SF1, double & SF2);

	void   Rotation							(Mat3_t Input, Mat3_t & Vectors, Mat3_t & VectorsT, Mat3_t & Values);

	Mat3_t abab									(Mat3_t const & A, Mat3_t const & B);

}; // namespace SPH

#include "Functions.cpp"

#endif // SPH_SPECIAL_FUNCTIONS_H
