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

#include "Functions.h"

namespace SPH {

inline double Kernel(size_t Dim, size_t KT, double r, double h)
{
	double C;
	double q = r/h;

	switch (KT)
	{
	case 0:
	    	// Qubic Spline
	    	Dim ==2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);

		if 	((q>=0.0)&&(q<1.0))	return C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
		else if (q<2.0)			return C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q));
		else				return 0.0;
		break;
	case 1:
	    	// Quadratic
	    	Dim ==2 ? C = 2.0/(h*h*M_PI) : C = 5.0/(4.0*h*h*h*M_PI);

	    	if	(q<2.0)			return C*(3.0/4.0-(3.0/4.0)*q+(3.0/16.0)*q*q);
		else				return 0.0;
	    	break;
	case 2:
	    	// Quintic
	    	Dim ==2 ? C = 7.0/(4.0*h*h*M_PI) : C = 7.0/(8.0*h*h*h*M_PI);

	    	if	(q<2.0)			return C*pow((1.0-q/2.0),4.0)*(2.0*q+1.0);
		else				return 0.0;
	    	break;
	case 3:
	    	// Gaussian with compact support
	    	Dim ==2 ? C = 1.0/(h*h*M_PI) : C = 1.0/(h*h*h*pow(M_PI,(3.0/2.0)));

	    	if (q<=2.0)			return C*exp(-q*q);
		else				return 0.0;
	    	break;
	case 4:
		// Quintic Spline
		Dim ==2 ? C = 7.0/(478.0*h*h*M_PI) : C = 1.0/(120.0*h*h*h*M_PI);

		if	((q>=0.0)&&(q<1.0))	return C*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0)+15.0*pow((1.0-q),5.0));
		else if (q<2.0)			return C*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0));
		else if (q<3.0)			return C*(pow((3.0-q),5.0));
		else				return 0.0;
		break;
	default:
	   	std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Qubic Spline" << std::endl;
		std::cout << "1 => Quadratic" << std::endl;
		std::cout << "2 => Quintic" << std::endl;
		std::cout << "3 => Gaussian with compact support of q<2" << std::endl;
		std::cout << "4 => Quintic Spline" << std::endl;
		abort();
		break;
	}
}

inline double GradKernel(size_t Dim, size_t KT, double r, double h)
{
	double C;
	double q = r/h;

	switch (KT)
	{
	case 0:
	    	// Qubic Spline
	    	Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);

		if 	(q==0.0)		return C/h*(-3.0+(9.0/2.0)*q);
		else if ((q>0.0)&&(q<1.0))	return C/r*(-3.0*q+(9.0/4.0)*q*q);
		else if (q<2.0)			return C/r*((-3.0/4.0)*(2.0-q)*(2.0-q));
		else				return 0.0;
		break;
	case 1:
	    	// Quadratic
	    	Dim ==2 ? C = 2.0/(h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*M_PI);

	    	if	(q<2.0)			return C/r*(-3.0/4.0+(3.0/8.0)*q);
		else				return 0.0;
	    	break;
	case 2:
	    	// Quintic
	    	Dim ==2 ? C = 7.0/(4.0*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*M_PI);

	    	if 	(q==0.0)		return C*-5.0/h*(pow((1.0-q/2.0),3.0)-3.0*q/2.0*pow((1.0-q/2.0),2.0));
	    	else if ((q>0.0)&&(q<2.0))	return C/r*-5.0*q*pow((1.0-q/2.0),3.0);
		else				return 0.0;
	    	break;
	case 3:
	    	// Gaussian with compact support
	    	Dim ==2 ? C = 1.0/(h*h*h*M_PI) : C = 1.0/(h*h*h*h*pow(M_PI,(3.0/2.0)));

	    	if 	(q==0.0)		return C*-2.0/h*(exp(-q*q)-2.0*q*q*exp(-q*q));
	    	else if ((q>0.0)&&(q<=2.0))	return C/r*-2.0*q*exp(-q*q);
		else				return 0.0;
	    	break;
	case 4:
		// Quintic Spline
		Dim ==2 ? C = 7.0/(478.0*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*M_PI);

		if	(q==0.0)		return C/h*(20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0)+300.0*pow((1.0-q),3.0));
		else if ((q>0.0)&&(q<1.0))	return C/r*(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0)-75.0*pow((1.0-q),4.0));
		else if (q<2.0)			return C/r*(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0));
		else if (q<3.0)			return C/r*(-5.0*pow((3.0-q),4.0));
		else				return 0.0;

		break;
	default:
	   	std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Qubic Spline" << std::endl;
		std::cout << "1 => Quadratic" << std::endl;
		std::cout << "2 => Quintic" << std::endl;
		std::cout << "3 => Gaussian with compact support of q<2" << std::endl;
		std::cout << "4 => Quintic Spline" << std::endl;
		abort();
		break;
	}
}

inline double LaplaceKernel(size_t Dim, size_t KT, double r, double h)
{
	double C;
	double q = r/h;

	switch (KT)
	{
	case 0:
		// Qubic Spline
	    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

		if	((q>=0.0)&&(q<1.0))	return C*(-3.0+(9.0/2.0)*q)  + C*(Dim-1.0)/q * (-3.0*q+(9.0/4.0)*q*q);
		else if (q<2.0) 		return C*((3.0/2.0)*(2.0-q)) + C*(Dim-1.0)/q * ((-3.0/4.0)*(2.0-q)*(2.0-q));
		else				return 0.0;
		break;
	case 1:
	    	// Quadratic
	    	Dim ==2 ? C = 2.0/(h*h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*h*M_PI);

	    	if	(q<2.0)			return C*(-3.0/8.0) + C*(Dim-1.0)/q * (-3.0/4.0+(3.0/8.0)*q);
		else				return 0.0;
	    	break;
	case 2:
	    	// Quintic
	    	Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

	    	if 	(q<2.0)			return C*pow((1.0-q/2.0),2.0)*(10.0*q-5.0) + C*(Dim-1.0)/q * -5.0*q*pow((1.0-q/2.0),3.0);
		else				return 0.0;
	    	break;
	case 3:
	    	// Gaussian with compact support
	    	Dim ==2 ? C = 1.0/(h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*pow(M_PI,(3.0/2.0)));

	    	if 	(q<=2.0)		return C*2.0*(2.0*q*q-1.0)*exp(-q*q) + C*(Dim-1.0)/q * -2.0*q*exp(-q*q);
		else				return 0.0;
	    	break;
	case 4:
		// Quintic Spline
		Dim ==2 ? C = 7.0/(478.0*h*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*h*M_PI);

		if	((q>=0.0)&&(q<1.0))	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2-q),3.0)+300.0*pow((1-q),3.0)) + C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0)-75.0*pow((1.0-q),4.0));
		else if (q<2.0)			return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2-q),3.0))                      + C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0));
		else if (q<3.0)			return C*(20.0*pow((3.0-q),3.0))                                           + C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0));
		else				return 0.0;
		break;
	default:
	   	std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Qubic Spline" << std::endl;
		std::cout << "1 => Quadratic" << std::endl;
		std::cout << "2 => Quintic" << std::endl;
		std::cout << "3 => Gaussian with compact support of q<2" << std::endl;
		std::cout << "4 => Quintic Spline" << std::endl;
		abort();
		break;
	}
}

inline double SecDerivativeKernel(size_t Dim, size_t KT, double r, double h)
{
	double C;
	double q = r/h;

	switch (KT)
	{
	case 0:
	    	// Qubic Spline
	    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

		if 	((q>=0.0)&&(q<1.0))	return C*(-3.0+(9.0/2.0)*q);
		else if (q<2.0)			return C*((3.0/2.0)*(2.0-q));
		else				return 0.0;
		break;
	case 1:
	    	// Quadratic
	    	Dim ==2 ? C = 2.0/(h*h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*h*M_PI);

	    	if	(q<2.0)			return C*(-3.0/8.0);
		else				return 0.0;
	   	break;
	case 2:
	    	// Quintic
	    	Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

	    	if	(q<2.0)			return C*pow((1.0-q/2.0),2.0)*(10.0*q-5.0);
	    	else				return 0.0;
	    	break;
	case 3:
	    	// Gaussian with compact support
	    	Dim ==2 ? C = 1.0/(h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*pow(M_PI,(3.0/2.0)));

	    	if	(q<=2.0)		return C*2.0*(2.0*q*q-1.0)*exp(-q*q);
		else				return 0.0;
	    	break;
	case 4:
		// Quintic Spline
		Dim ==2 ? C = 7.0/(478.0*h*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*h*M_PI);

		if	((q>=0.0)&&(q<1.0))	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0)+300.0*pow((1.0-q),3.0));
		else if (q<2.0)			return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0));
		else if (q<3.0)			return C*(20.0*pow((3.0-q),3.0));
		else				return 0.0;
		break;
	default:
	   	std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Qubic Spline" << std::endl;
		std::cout << "1 => Quadratic" << std::endl;
		std::cout << "2 => Quintic" << std::endl;
		std::cout << "3 => Gaussian with compact support of q<2" << std::endl;
		std::cout << "4 => Quintic Spline" << std::endl;
		abort();
		break;
	}
}

inline double EOS(size_t EQ, double Cs0, double P00, double Density, double Density0)
{
	switch (EQ)
	{
	case 0:
		return P00+(Cs0*Cs0)*(Density-Density0);
		break;
	case 1:
		return P00+(Density0*Cs0*Cs0/7.0)*(pow(Density/Density0,7.0)-1);
		break;
	case 2:
		return (Cs0*Cs0)*Density;
		break;
	default:
		std::cout << "Please correct Pressure Equation No and run again" << std::endl;
		std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
		std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
		std::cout << "2 => (Cs*Cs)*Density" << std::endl;
		abort();
		break;
	}
}

inline double SoundSpeed(size_t EQ, double Cs0, double Density, double Density0)
{
	switch (EQ)
	{
	case 0:
		return Cs0;
		break;
	case 1:
		return sqrt((Cs0*Cs0)*pow(Density/Density0,6.0));
		break;
	case 2:
		return Cs0;
		break;
	default:
		std::cout << "Please correct Pressure Equation No and run again" << std::endl;
		std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
		std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
		std::cout << "2 => (Cs*Cs)*Density" << std::endl;
		abort();
		break;
	}
}

inline double DensitySolid(size_t EQ, double Cs0, double P00, double Pressure, double Density0)
{
	switch (EQ)
	{
	case 0:
		return (Pressure-P00)/(Cs0*Cs0) + Density0;
		break;
	case 1:
		return pow( ((Pressure-P00)*(7.0/(Density0*Cs0*Cs0))+1) , 1.0/7.0 ) * Density0;
		break;
	case 2:
		return Pressure/(Cs0*Cs0);
		break;
	default:
		std::cout << "Please correct Pressure Equation No and run again" << std::endl;
		std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
		std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
		std::cout << "2 => (Cs*Cs)*Density" << std::endl;
		abort();
		break;
	}
}

inline void Seepage(size_t ST, double k, double k2, double mu,  double rho, double& SF1, double& SF2)
{
	switch (ST)
	{
	case 0:
		// Darcy
		SF1 = mu/k;
		SF2 = 0.0;
		break;
	case 1:
		// Darcy_Kozeny–Carman EQ
		SF1 = mu/k;
		SF2 = 0.0;
		break;
	case 2:
		// Ergun
		SF1 = mu/k;
		SF2 = k2*rho;
		break;
	case 3:
		// Den Adel
		SF1 = mu/k;
		SF2 = k2*rho;
		break;
	default:
		std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Darcy's Law" << std::endl;
		std::cout << "1 => Darcy's Law & Kozeny–Carman Eq" << std::endl;
		std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
		std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
		abort();
		break;
	}
}

inline void Rotation (Mat3_t Input, Mat3_t & Vectors, Mat3_t & VectorsT, Mat3_t & Values)
{
	Vec3_t Val,V0,V1,V2;
	Eig(Input,Val,V0,V1,V2,true,false);

	Mat3_t V;
	Vectors	=	V0(0), V1(0), V2(0),
			V0(1), V1(1), V2(1),
			V0(2), V1(2), V2(2);

	Trans(Vectors,VectorsT);
	Values	= 	Val(0), 0.0   , 0.0   ,
			0.0   , Val(1), 0.0   ,
			0.0   , 0.0   , Val(2);
}

inline Mat3_t abab (const Mat3_t & A, const Mat3_t & B)
{
		Mat3_t M;
		M(0,0)=A(0,0)*B(0,0);  M(0,1)=A(0,1)*B(0,1);  M(0,2)=A(0,2)*B(0,2);
		M(1,0)=A(1,0)*B(1,0);  M(1,1)=A(1,1)*B(1,1);  M(1,2)=A(1,2)*B(1,2);
		M(2,0)=A(2,0)*B(2,0);  M(2,1)=A(2,1)*B(2,1);  M(2,2)=A(2,2)*B(2,2);
		return M;
}

}; // namespace SPH
