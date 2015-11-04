/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo Torres  *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_SPH_PARTICLE_H
#define MECHSYS_SPH_PARTICLE_H

// Std lib
#include <iostream>
#include <cstring>
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
    bool   	IsFree;			///< Check the particle if it is free to move or not
    bool   	NoSlip;			///< No-Slip BC
    Vec3_t	FreeSlip;
    int    	ID;				///< an Integer value to identify the particle set
    int    	Material;		///< an Integer value to identify the particle material type
    						///< 1 = Fluid, 2 = Solid

    Vec3_t  x;				///< Position of the particle n
    Vec3_t  vb;				///< Velocity of the particle n-2
    Vec3_t  v;				///< Velocity of the particle n,
    Vec3_t	VXSPH;			///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
    Vec3_t  a;				///< Acceleration of the particle
    Vec3_t  Vis;			///< Acceleration of the particle due to viscosity effect

    double	ZWab;			///< Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
    double	SumDen;			///< Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)

    double 	Pressure;		///< Pressure of the particle n

    double	Density;		///< Density of the particle n
    double 	Densityb;		///< Density of the particle n-2
    double 	RefDensity;		///< Reference Density of Particle
    double 	dDensity;		///< Rate of density change in time based on state equations

    double 	Mass;			///< Mass of the particle

    Mat3_t  StrainRate;		///< Global shear Strain rate tensor
    Mat3_t  Rotation;		///< Global rotation tensor

    double  ShearRate;		///< Global shear rate

    Mat3_t  ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n
    Mat3_t  ShearStressb;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-2
    Mat3_t  Sigma;			///< Cauchy stress tensor (Total Stress) n
    Mat3_t  Sigmab;			///< Cauchy stress tensor (Total Stress) n-2
    Mat3_t  TIR;			///< Tensile Instability stress tensor R
    Mat3_t  Strain;			///< Total Strain n
    Mat3_t  Strainb;		///< Total Strain n-2

    double 	Mu;				///< Dynamic viscosity coefficient of the fluid particle
    double 	MuRef;			///< Reference Dynamic viscosity coefficient
    double 	T0;		  		///< Yield stress for Bingham fluids
    double 	m;		  		///< Normalization value for Bingham fluids

    double 	G;				///< Shear modulus
    double 	K;				///< Bulk modulus

    size_t	Fail;			///< Failure criteria
    double	c;				///< Cohesion
    double	phi;			///< Friction angel
    double	psi;			///< Dilation angel
    double	Sigmay;			///< Tensile yield stress


    double 	h;				///< Smoothing length of the particle
    double 	hr;				///< Reference smoothing length of the particle

    int    	LL;				///< Linked-List variable to show the next particle in the list of a cell
    int    	CC[3];			///< Current cell No for the particle (linked-list)

    int		ct;				///< Correction step for the Verlet Algorithm and Shepard filter

    omp_lock_t my_lock;		///< Open MP lock

    double SumKernel;		///<Summation of the kernel value for neighbour particles

    // Constructor
    Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed=false);

    // Methods
    void Move			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin,
    						bool ShepardFilter, size_t PresEq, double Cs, double P0);	///< Update the important quantities of a particle
    void translate		(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
    void Mat1			(double dt, bool ShepardFilter, size_t PresEq, double Cs, double P0);
    void Mat2			(double dt, size_t PresEq, double Cs, double P0);
    void Mat3			(double dt);
};

inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0,bool Fixed)
{
	ct = 0;
	a = 0.0;
    x = x0;
    vb = v = v0;
    Densityb = Density = Density0;
    RefDensity = Density0;
    Mass = Mass0;
    IsFree = !Fixed;
    hr = h0;
    h = hr;
    Pressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    LL=0;
    omp_init_lock(&my_lock);
    VXSPH = 0.0;
    ZWab = 0.0;
    SumDen = 0.0;
    dDensity=0.0;
    ShearRate = 0.0;
    StrainRate = 0.0;
    Rotation = 0.0;
    MuRef = Mu = 0.0;
    T0 = 0.0;
    m = 300.0;
    SumKernel = 0.0;
    G = 0.0;
    Material = 0;
    Fail = 0;
    c = 0.0;
    phi = 0.0;
    psi = 0.0;
    Sigmay = 0.0;
    K = 0.0;
    NoSlip = false;
    FreeSlip = 0.0;
    set_to_zero(ShearStress);
    set_to_zero(ShearStressb);
    set_to_zero(TIR);

}

inline void Particle::Move (double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, bool ShepardFilter, size_t PresEq, double Cs, double P0)
{
	// Evolve position
	x = x + dt*(v+VXSPH) + 0.5*dt*dt*a;

	// Evolve velocity
	if (ct==30)
	{
		Vec3_t temp;
		temp = v;
		v = v + dt*a;
		vb = temp;
		ct = -1;
	}
	else
	{
		Vec3_t temp;
		temp = v;
		v = vb + 2*dt*a;
		vb = temp;
	}
	ct++;

	switch (Material)
    {case 1:
    	Mat1(dt,ShepardFilter,PresEq,Cs,P0);
		break;
    case 2:
    	Mat2(dt,PresEq,Cs,P0);
    	break;
    case 3:
    	Mat3(dt);
    	break;
   default:
	   	std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
	    abort();
	    break;
    }

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

inline void Particle::Mat1(double dt, bool ShepardFilter, size_t PresEq, double Cs, double P0)
{
	// Evolve density
	if (ShepardFilter && ct==30)
	{
		if (!isnan(SumDen/ZWab))
		{
			// Shepard filter
			Density = SumDen/ZWab;
			Densityb = Density;
		}
		else
			Densityb = Density;
	}
	else
	{
		double dens = Density;
		Density = Density + dt*dDensity;
		Densityb = dens;
	}

	Pressure = PressureEq(PresEq, Cs, P0,Density, RefDensity);

	ShearRate = sqrt(0.5*(StrainRate(0,0)*StrainRate(0,0) + 2.0*StrainRate(0,1)*StrainRate(1,0) +
			2.0*StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,1)*StrainRate(1,1) +
			2.0*StrainRate(1,2)*StrainRate(2,1) + StrainRate(2,2)*StrainRate(2,2)));

	// Bingham viscosity calculation
	if (T0>0.0)
	{
		if (ShearRate !=0.0)
			Mu = MuRef + T0*(1-exp(-m*ShearRate))/ShearRate;
		else
			Mu = MuRef + T0*m;
	}

}

inline void Particle::Mat2(double dt, size_t PresEq, double Cs, double P0)
{
	// Evolve density
	double dens = Density;
	Density = Density + dt*dDensity;
	Densityb = dens;

	Pressure = PressureEq(PresEq, Cs, P0,Density, RefDensity);

	Mat3_t RotationRateT, Stress,SRT,RS;

	// Jaumann rate terms
	Trans(Rotation,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(Rotation,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	Stress			= ShearStress;
	ShearStress		= 2.0*dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressb;
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

	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}

inline void Particle::Mat3(double dt)
{
	Mat3_t RotationRateT, Stress, SRT,RS;
	double I1,J2,alpha,k,I1strain;

	Stress	= Strain;
	Strain	= 2.0*dt*StrainRate + Strainb;
	Strainb	= Stress;

	// Jaumann rate terms
	Trans(Rotation,RotationRateT);
	Mult(Sigma,RotationRateT,SRT);
	Mult(Rotation,Sigma,RS);

	// Volumetric strain
	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2);

	// Elastic prediction step (Sigma_e n+1)
	Stress	= Sigma;
	Sigma	=	2.0*dt*(
							I1strain*K*OrthoSys::I
						+	2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I)
						+	SRT+RS
						) + Sigmab;
	Sigmab	= Stress;

	if (Fail>=2)
	{
		// Drucker-Prager failure criterion for plane strain
		I1		= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
		alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
		k		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));

		// Bring back stress to the apex of the failure criteria
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
		if ((sqrt(J2)+alpha*I1-k)>=0.0 && sqrt(J2)>0.0)
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
					Plastic	= 3.0*alpha*K*OrthoSys::I + G/sqrt(J2)*ShearStress;
					break;
				case 3:
					dLanda	= 1.0/(9.0*alpha*K*3.0*sin(psi)+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
					Plastic	= 3.0*3.0*sin(psi)*K*OrthoSys::I + G/sqrt(J2)*ShearStress;
					break;
				default:
					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "2 => Associated flow rule" << std::endl;
					std::cout << "3 => non-associated flow rule" << std::endl;
					abort();
					break;
				}
				// Apply the plastic term
			    if (dLanda>0.0) Sigma = Sigma -	2.0*dt*(fabs(dLanda)*Plastic);
			}
			//Scale back
			I1			= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
			ShearStress	= Sigma - 1.0/3.0* I1 *OrthoSys::I;
			J2			= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
							2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
							2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));

			if ((sqrt(J2)+alpha*I1-k)>0.0 && sqrt(J2)>0.0) Sigma = I1/3.0*OrthoSys::I + (k-alpha*I1)/sqrt(J2) * ShearStress;
		}
	}
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
