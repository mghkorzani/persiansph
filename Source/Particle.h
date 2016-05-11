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

#ifndef MECHSYS_SPH_PARTICLE_H
#define MECHSYS_SPH_PARTICLE_H

// Std lib
#include <iostream>
#include <cmath>
#include <omp.h>

// MechSys
#include "../External/matvec.h"
#include "Functions.h"

namespace SPH {
Mat3_t abab (const Mat3_t & A, const Mat3_t & B)
{
    Mat3_t M;
    M(0,0)=A(0,0)*B(0,0);  M(0,1)=A(0,1)*B(0,1);  M(0,2)=A(0,2)*B(0,2);
    M(1,0)=A(1,0)*B(1,0);  M(1,1)=A(1,1)*B(1,1);  M(1,2)=A(1,2)*B(1,2);
    M(2,0)=A(2,0)*B(2,0);  M(2,1)=A(2,1)*B(2,1);  M(2,2)=A(2,2)*B(2,2);
    return M;
}

class Particle
{
public:
    // Data
    bool   	Shepard;	///< Shepard Filter for the density
    size_t	ShepardCounter;	///< Count number of contributing particles
    size_t	ShepardStep;
    size_t	ShepardNeighbourNo;
    bool   	IsFree;		///< Check the particle if it is free to move or not
    bool   	InOut;		///< Check the particle if it is in-flow or out-flow or not
    bool   	IsSat;		///< Check the particle if it is Saturated or not
    bool   	SatCheck;	///< Check the particle Saturation at each time step
    bool   	NoSlip;		///< No-Slip BC
    int    	ID;				///< an Integer value to identify the particle set
    int    	Material;	///< an Integer value to identify the particle material type: 1 = Fluid, 2 = Solid, 3 = Soil

    Vec3_t  x;			///< Position of the particle n
    Vec3_t  vb;			///< Velocity of the particle n-1 (Modified Verlet)
    Vec3_t  va;			///< Velocity of the particle n+1/2 (Leapfrog)
    Vec3_t  v;			///< Velocity of the particle n+1,
    Vec3_t	VXSPH;		///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
    Vec3_t  a;			///< Acceleration of the particle n
    Vec3_t  Vis;		///< Acceleration of the particle due to viscosity effect n

    double	ZWab;		///< Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
    double	SumDen;		///< Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)

    double 	Pressure;	///< Pressure of the particle n+1
    double	Density;	///< Density of the particle n+1
    double 	Densitya;	///< Density of the particle n+1/2 (Leapfrog)
    double 	Densityb;	///< Density of the particle n-1 (Modified Verlet)
    double 	dDensity;	///< Rate of density change in time based on state equations n

    double	V;			///< Volume of a particle
    double	RhoF;		///< Density of water or any other fluids

    double 	RefDensity;	///< Reference Density of Particle
    double 	Mass;		///< Mass of the particle

    Mat3_t  StrainRate;		///< Global shear Strain rate tensor n
    Mat3_t  RotationRate;		///< Global rotation tensor n
    double  ShearRate;		///< Global shear rate for fluids

    Mat3_t  ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1
    Mat3_t  ShearStressa;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1/2 (Leapfrog)
    Mat3_t  ShearStressb;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-1 (Modified Verlet)
    Mat3_t  Sigma;		///< Cauchy stress tensor (Total Stress) n+1
    Mat3_t  Sigmaa;		///< Cauchy stress tensor (Total Stress) n+1/2 (Leapfrog)
    Mat3_t  Sigmab;		///< Cauchy stress tensor (Total Stress) n-1 (Modified Verlet)

    Mat3_t  TIR;		///< Tensile Instability stress tensor R

    Mat3_t  Strain;		///< Total Strain n+1
    Mat3_t  Straina;		///< Total Strain n+1/2 (Leapfrog)
    Mat3_t  Strainb;		///< Total Strain n-1 (Modified Verlet)

    double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
    double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
    double 	Mu;		///< Dynamic viscosity coefficient of the fluid particle
    double 	MuRef;		///< Reference Dynamic viscosity coefficient
    double 	T0;		///< Yield stress for Bingham fluids
    double 	m;		///< Normalization value for Bingham fluids
    size_t	VisM;		///< Non-Newtonian viscosity method
	
    double 	G;		///< Shear modulus
    double 	K;		///< Bulk modulus

    double	TI;		///< Tensile instability factor
    double	TIn;		///< Tensile instability power

    size_t	Fail;		///< Failure criteria
    double	c;				///< Cohesion
    double	phi;		///< Friction angel
    double	psi;		///< Dilation angel
    double	Sigmay;		///< Tensile yield stress
    double	n;		///< Prosity
    double	k;		///< Permeability
    double	d;		///< effective particle size

    double 	h;		///< Smoothing length of the particle

    int    	LL;		///< Linked-List variable to show the next particle in the list of a cell
    int    	CC[3];		///< Current cell No for the particle (linked-list)

    int		ct;		///< Correction step for the Modified Verlet Algorithm and Shepard filter

    omp_lock_t my_lock;		///< Open MP lock

    double	SumKernel;	///< Summation of the kernel value for neighbour particles
    bool	FirstStep;	///< to initialize the integration scheme
    size_t	PresEq;		///< Selecting variable to choose an equation of state
    double	Cs;		///< Speed of sound
    double	P0;		///< background pressure for equation of state

    // Constructor
    Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed=false);

    // Methods
    void Move			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin,
    						size_t Scheme, Mat3_t I);	///< Update the important quantities of a particle
    void Move_MVerlet	(Mat3_t I, double dt);	///< Update the important quantities of a particle
    void Move_Leapfrog	(Mat3_t I, double dt);	///< Update the important quantities of a particle
    void translate		(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
    void Mat1			(double dt);
    void Mat2MVerlet	(double dt);
    void Mat3MVerlet	(Mat3_t I, double dt);
    void Mat2Leapfrog	(double dt);
    void Mat3Leapfrog	(Mat3_t I, double dt);
};

inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0,bool Fixed)
{
	ct = 0;
	a = 0.0;
    x = x0;
    n = 0.0;
    k = 0.0;

    Cs		= 0.0;
    P0		= 0.0;
    PresEq	= 0;
    Alpha	= 0.0;
    Beta	= 0.0;

    va = 0.0;
    vb = 0.0;
    v = v0;
    VXSPH = 0.0;
    TI		= 0.0;
    TIn		= 4.0;

    Densitya = 0.0;
    Densityb = 0.0;
    Density = Density0;
    RefDensity = Density0;

    Mass = Mass0;
    IsFree = !Fixed;
    h = h0;
    Pressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    LL=0;
    ZWab = 0.0;
    SumDen = 0.0;
    dDensity=0.0;
    ShearRate = 0.0;
    MuRef = Mu = 0.0;
    VisM = 0;
    T0 = 0.0;
    m = 300.0;
    SumKernel = 0.0;
    G = 0.0;
    K = 0.0;
    Material = 0;
    Fail = 0;
    c = 0.0;
    phi = 0.0;
    psi = 0.0;
    d =0.0;
    Sigmay = 0.0;
    NoSlip = false;
    Shepard = false;
    InOut = false;
    FirstStep = true;
    V = Mass/RefDensity;
    RhoF = 0.0;
    IsSat = false;
    SatCheck = false;
    ShepardNeighbourNo = 0;
    ShepardStep = 40;
    ShepardCounter = 0;


    set_to_zero(Strainb);
    set_to_zero(Strain);
    set_to_zero(Sigmab);
    set_to_zero(Sigma);
    set_to_zero(Sigmaa);
    set_to_zero(ShearStress);
    set_to_zero(ShearStressb);
    set_to_zero(TIR);
    set_to_zero(StrainRate);
    set_to_zero(RotationRate);
    omp_init_lock(&my_lock);

}

inline void Particle::Move(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, size_t Scheme, Mat3_t I)
{
	if (Scheme == 0)
		Move_MVerlet(I, dt);
	else
		Move_Leapfrog(I, dt);


	//Periodic BC particle position update
	if (Domainsize(0)>0.0)
	{
		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1)>0.0)
	{
		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2)>0.0)
	{
		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}

}

inline void Particle::Mat1(double dt)
{
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	ShearRate = sqrt(0.5*(StrainRate(0,0)*StrainRate(0,0) + 2.0*StrainRate(0,1)*StrainRate(1,0) +
			2.0*StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,1)*StrainRate(1,1) +
			2.0*StrainRate(1,2)*StrainRate(2,1) + StrainRate(2,2)*StrainRate(2,2)));

	// Bingham viscosity calculation
	if (T0>0.0)
	{
		switch (VisM)
		{
			case 0:
			// Bingham
				if (ShearRate !=0.0)
					Mu = MuRef + T0*(1-exp(-m*ShearRate))/ShearRate;
				else
					Mu = MuRef + T0*m;
				break;
			case 1:
			// Cross
				Mu = (1000.0*MuRef + MuRef^2*1000.0/T0*ShearRate)/(1+1000.0*MuRef/T0*ShearRate);
				break;
			default:
				std::cout << "Non-Newtonian Viscosity Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => Bingham" << std::endl;
				std::cout << "1 => Cross" << std::endl;
				abort();
				break;
		}
	}
}


inline void Particle::Move_MVerlet (Mat3_t I, double dt)
{
	if (FirstStep)
	{
		ct = 30;
		FirstStep = false;
	}

	x += dt*(v+VXSPH) + 0.5*dt*dt*a;

	if (ct == 30)
	{
		if (Shepard && ShepardCounter == ShepardStep && Material == 1)
		{
			if (ShepardNeighbourNo>=3)
			{
				Densityb	= SumDen/ZWab;
//				Densityb	= Density;
				Density		= SumDen/ZWab;
			}
			else
			{
				Densityb	= Density;
				Density		+=dt*dDensity;
			}
			ShepardNeighbourNo = 0;
		}
		else
		{
			Densityb		= Density;
			Density			+=dt*dDensity;
		}

		vb	= v;
		v	+=dt*a;
	}
	else
	{
		if (Shepard && ShepardCounter == ShepardStep && Material == 1)
		{
			if (ShepardNeighbourNo>=3)
			{
				Densityb	= SumDen/ZWab;
//				Densityb	= Density;
				Density		= SumDen/ZWab;
			}
			else
			{
				double dens	= Density;
				Density		= Densityb + 2.0*dt*dDensity;
				Densityb	= dens;
			}
			ShepardNeighbourNo = 0;
		}
		else
		{
			double dens	= Density;
			Density		= Densityb + 2.0*dt*dDensity;
			Densityb	= dens;
		}

		Vec3_t temp;
		temp	= v;
		v		= vb + 2*dt*a;
		vb		= temp;
	}

	switch (Material)
    {case 1:
    	Mat1(dt);
		break;
    case 2:
    	Mat2MVerlet(dt);
    	break;
    case 3:
    	Mat3MVerlet(I,dt);
    	break;
   default:
	   	std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
	    abort();
	    break;
    }
	if (ct == 30) ct = 0; else ct++;
	if (ShepardCounter == ShepardStep) ShepardCounter = 0; else ShepardCounter++;
}

inline void Particle::Mat2MVerlet(double dt)
{
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT, Stress,SRT,RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	Stress			= ShearStress;
	if (ct == 30)
		ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
	else
		ShearStress	= 2.0*dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressb;
	ShearStressb	= Stress;

	if (Fail == 1)
	{
		double J2	= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));
		//Scale back
		ShearStress	= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStress;
	}

	Sigma			= -Pressure * OrthoSys::I + ShearStress;

	Stress	= Strain;
	if (ct == 30)
		Strain	= dt*StrainRate + Strain;
	else
		Strain	= 2.0*dt*StrainRate + Strainb;
	Strainb	= Stress;


	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}

inline void Particle::Mat3MVerlet(Mat3_t I, double dt)
{
	Mat3_t RotationRateT, Stress, SRT,RS;
	double I1,J2,alpha,k,I1strain;

	// Jaumann rate terms
	Trans(RotationRate,RotationRateT);
	Mult(Sigma,RotationRateT,SRT);
	Mult(RotationRate,Sigma,RS);

	// Volumetric strain
	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2);

	// Elastic prediction step (Sigma_e n+1)
	Stress	= Sigma;
	if (ct == 30)
		Sigma	= dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigma;
	else
		Sigma	= 2.0*dt*( I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigmab;
	Sigmab	= Stress;

	if (Fail>1)
	{
		// Drucker-Prager failure criterion for plane strain
		alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
		k		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));

		// Bring back stress to the apex of the failure criteria
		I1		= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
		if ((k-alpha*I1)<0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio =0.0; else Ratio = k/alpha;
			Sigma(0,0) -= 1.0/3.0*(I1-Ratio);
			Sigma(1,1) -= 1.0/3.0*(I1-Ratio);
			Sigma(2,2) -= 1.0/3.0*(I1-Ratio);
			I1 			= Ratio;
		}

		// Shear stress based on the elastic assumption (S_e n+1)
		ShearStress = Sigma - 1.0/3.0* I1 *OrthoSys::I;
		J2 			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));


		// Check the elastic prediction step by the failure criteria
		if ((sqrt(J2)+alpha*I1-k)>0.0)
		{
			// Shear stress based on the existing stress (S n)
			ShearStress = Stress - 1.0/3.0*(Stress(0,0)+Stress(1,1)+Stress(2,2))*OrthoSys::I;
			J2 			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
							2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
							2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));

			if (sqrt(J2)>0.0)
			{
				Mat3_t temp, Plastic;
				double sum,dLanda;

				// calculating the plastic term based on the existing shear stress and strain rate
				temp	= abab(ShearStress,StrainRate);
				sum		= temp(0,0)+temp(0,1)+temp(0,2)+temp(1,0)+temp(1,1)+temp(1,2)+temp(2,0)+temp(2,1)+temp(2,2);
				switch (Fail)
				{
				case 2:
					dLanda	= 1.0/(9.0*alpha*alpha*K+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
					Plastic	= 3.0*alpha*K*I + G/sqrt(J2)*ShearStress;
					break;
				case 3:
					dLanda	= 1.0/(9.0*alpha*K*3.0*sin(psi)+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
					Plastic	= 3.0*3.0*sin(psi)*K*I + G/sqrt(J2)*ShearStress;
					break;
				default:
					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "2 => Associated flow rule" << std::endl;
					std::cout << "3 => non-associated flow rule" << std::endl;
					abort();
					break;
				}
				// Apply the plastic term
				if (ct == 30)
					Sigma = Sigma -	dt*(dLanda*Plastic);
				else
					Sigma = Sigma -	2.0*dt*(dLanda*Plastic);
			}

			//Scale back
			I1			= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
			if ((k-alpha*I1)<0.0)
			{
				double Ratio;
				if (alpha == 0.0) Ratio =0.0; else Ratio = k/alpha;
				Sigma(0,0) -= 1.0/3.0*(I1-Ratio);
				Sigma(1,1) -= 1.0/3.0*(I1-Ratio);
				Sigma(2,2) -= 1.0/3.0*(I1-Ratio);
				I1 			= Ratio;
			}
			ShearStress	= Sigma - 1.0/3.0* I1 *OrthoSys::I;
			J2			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
							2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
							2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));

			if ((sqrt(J2)+alpha*I1-k)>0.0 && sqrt(J2)>0.0) Sigma = I1/3.0*OrthoSys::I + (k-alpha*I1)/sqrt(J2) * ShearStress;
		}
	}

	Stress	= Strain;
	if (ct == 30)
		Strain	= dt*StrainRate + Strain;
	else
		Strain	= 2.0*dt*StrainRate + Strainb;
	Strainb	= Stress;

}

inline void Particle::Move_Leapfrog(Mat3_t I, double dt)
{
	if (FirstStep)
	{
		Densitya = Density - dt/2.0*dDensity;
		va = v - dt/2.0*a;
	}
	Densityb = Densitya;
	Densitya += dt*dDensity;
	Density = (Densitya+Densityb)/2.0;
	vb = va;
	va += dt*a;
	v = (va + vb)/2.0;
	x += dt*va;

	switch (Material)
    {case 1:
    	Mat1(dt);
		break;
    case 2:
    	Mat2Leapfrog(dt);
    	break;
    case 3:
    	Mat3Leapfrog(I,dt);
    	break;
   default:
	   	std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
	    abort();
	    break;
    }
	if (FirstStep) FirstStep = false;

}

inline void Particle::Mat2Leapfrog(double dt)
{
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT,SRT,RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	if (FirstStep)
		ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

	ShearStressb	= ShearStressa;
	ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;

	if (Fail == 1)
	{
		double J2	= 0.5*(ShearStressa(0,0)*ShearStressa(0,0) + 2.0*ShearStressa(0,1)*ShearStressa(1,0) +
						2.0*ShearStressa(0,2)*ShearStressa(2,0) + ShearStressa(1,1)*ShearStressa(1,1) +
						2.0*ShearStressa(1,2)*ShearStressa(2,1) + ShearStressa(2,2)*ShearStressa(2,2));
		//Scale back
		ShearStressa= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;
	}
	ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);

	Sigma = -Pressure * OrthoSys::I + ShearStress;

	if (FirstStep)
		Straina	= -dt/2.0*StrainRate + Strain;
	Strainb	= Straina;
	Straina	= dt*StrainRate + Straina;
	Strain	= 1.0/2.0*(Straina+Strainb);


	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}

inline void Particle::Mat3Leapfrog(Mat3_t I, double dt)
{
	Mat3_t RotationRateT, Stress, SRT,RS;
	double I1,J2,alpha,k,I1strain;

	// Jaumann rate terms
	Trans(RotationRate,RotationRateT);
	Mult(Sigma,RotationRateT,SRT);
	Mult(RotationRate,Sigma,RS);

	// Volumetric strain
	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2);

	// Elastic prediction step (Sigma_e n+1)
	if (FirstStep)
		Sigmaa	= -dt/2.0*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigma;

	Sigmab	= Sigmaa;
	Sigmaa	= dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigmaa;

	if (Fail>1)
	{
		// Drucker-Prager failure criterion for plane strain
		alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
		k		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));

		// Bring back stress to the apex of the failure criteria
		I1		= Sigmaa(0,0) + Sigmaa(1,1) + Sigmaa(2,2);
		if ((k-alpha*I1)<0.0)
		{
			double Ratio;
			if (alpha == 0.0) Ratio =0.0; else Ratio = k/alpha;
			Sigmaa(0,0) -= 1.0/3.0*(I1-Ratio);
			Sigmaa(1,1) -= 1.0/3.0*(I1-Ratio);
			Sigmaa(2,2) -= 1.0/3.0*(I1-Ratio);
			I1 			= Ratio;
		}

		// Shear stress based on the elastic assumption (S_e n+1)
		ShearStress = Sigmaa - 1.0/3.0* I1 *OrthoSys::I;
		J2 			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));


		// Check the elastic prediction step by the failure criteria
		if ((sqrt(J2)+alpha*I1-k)>0.0)
		{
			// Shear stress based on the existing stress (S n)
			ShearStress = Sigma - 1.0/3.0*(Sigma(0,0)+Sigma(1,1)+Sigma(2,2))*OrthoSys::I;
			J2 			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
							2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
							2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));

			if (sqrt(J2)>0.0)
			{
				Mat3_t temp, Plastic;
				double sum,dLanda;

				// calculating the plastic term based on the existing shear stress and strain rate
				temp	= abab(ShearStress,StrainRate);
				sum		= temp(0,0)+temp(0,1)+temp(0,2)+temp(1,0)+temp(1,1)+temp(1,2)+temp(2,0)+temp(2,1)+temp(2,2);
				switch (Fail)
				{
				case 2:
					dLanda	= 1.0/(9.0*alpha*alpha*K+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
					Plastic	= 3.0*alpha*K*I + G/sqrt(J2)*ShearStress;
					break;
				case 3:
					dLanda	= 1.0/(9.0*alpha*K*3.0*sin(psi)+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
					Plastic	= 3.0*3.0*sin(psi)*K*I + G/sqrt(J2)*ShearStress;
					break;
				default:
					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "2 => Associated flow rule" << std::endl;
					std::cout << "3 => non-associated flow rule" << std::endl;
					abort();
					break;
				}
				Sigmaa = Sigmaa - dt*(dLanda*Plastic);
			}

			I1	= Sigmaa(0,0) + Sigmaa(1,1) + Sigmaa(2,2);
			if ((k-alpha*I1)<0.0)
			{
				double Ratio;
				if (alpha == 0.0) Ratio =0.0; else Ratio = k/alpha;
				Sigmaa(0,0) -= 1.0/3.0*(I1-Ratio);
				Sigmaa(1,1) -= 1.0/3.0*(I1-Ratio);
				Sigmaa(2,2) -= 1.0/3.0*(I1-Ratio);
				I1 			= Ratio;
			}
			ShearStress	= Sigmaa - 1.0/3.0* I1 *OrthoSys::I;
			J2			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
							2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
							2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));
			if ((sqrt(J2)+alpha*I1-k)>0.0 && sqrt(J2)>0.0) Sigmaa = I1/3.0*OrthoSys::I + (k-alpha*I1)/sqrt(J2) * ShearStress;
		}
	}
	Sigma = 1.0/2.0*(Sigmaa+Sigmab);

	if (FirstStep)
		Straina	= -dt/2.0*StrainRate + Strain;
	Strainb	= Straina;
	Straina	= dt*StrainRate + Straina;
	Strain	= 1.0/2.0*(Straina+Strainb);
}

inline void Particle::translate(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
{
	x = x + dt*v + 0.5*dt*dt*a;

	// Evolve velocity
	Vec3_t temp;
	temp = v;
	v = vb + 2*dt*a;
	vb = temp;

	//Periodic BC particle position update
	if (Domainsize(0)>0.0)
	{
		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1)>0.0)
	{
		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2)>0.0)
	{
		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}
}

}; // namespace SPH

#endif // MECHSYS_SPH_PARTICLE_H
