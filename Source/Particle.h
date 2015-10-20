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

namespace SPH {

class Particle
{
public:
    // Data
    bool   	IsFree;			///< Check the particle if it is free to move or not
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

    double 	Pressure;		///< Pressure at position of the particle

    double	Density;		///< Density at the position of the particle n
    double 	Densityb;		///< Density at the position of the particle n-2
    double 	RefDensity;		///< Reference Density of Particle
    double 	dDensity;		///< Rate of density change in time based on state equations

    double 	Mass;			///< Mass of the particle

    Mat3_t  StrainRate;		///< Global shear Strain rate tensor
    Mat3_t  Rotation;		///< Global rotation tensor

    double  ShearRate;		///< Global shear rate

    Mat3_t  ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor)
    double  NormalStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor)
    Mat3_t  Sigma;			///< Cauchy stress tensor

    double 	Mu;				///< Dynamic viscosity coefficient of the fluid particle
    double 	MuRef;			///< Reference Dynamic viscosity coefficient
    double 	T0;		  		///< Yield stress for Bingham fluids
    double 	m;		  		///< Normalization value for Bingham fluids

    double 	G;				///< Shear modulus
    double 	K;				///< Shear modulus

    size_t	Fail;			///< Failure criteria
    double	c;				///< Cohesion
    double	phi;			///< Friction angel
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
    void Move			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, bool ShepardFilter, Mat3_t I);	///< Update the important quantities of a particle
    void translate		(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
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
    Sigmay = 0.0;
    NormalStress = 0.0;
    K = 0.0;
}

inline void Particle::Move (double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, bool ShepardFilter, Mat3_t I)
{
	if (ct < 30)
	{
		if (IsFree)
		{
			// Evolve position
			x = x + dt*(v+VXSPH) + 0.5*dt*dt*a;

			// Evolve velocity
			Vec3_t temp;
			temp = v;
			v = vb + 2*dt*a;
			vb = temp;

			// Evolve density
			double dens = Density;
			Density = Densityb + 2*dt*dDensity;
			Densityb = dens;

			// Evolve shear stress
			if (Material == 2)
			{
				Mat3_t RotationRateT;
				Mat3_t SRT,RS;
				Trans(Rotation,RotationRateT);
				Mult(ShearStress,RotationRateT,SRT);
				Mult(Rotation,ShearStress,RS);
				ShearStress = dt*(2.0*G*(StrainRate-1.0/(2.0+I(2,2))*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*I)+SRT+RS) + ShearStress;
				NormalStress = dt*(1.0/(2.0+I(2,2))*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*K) + NormalStress;
			}
		}
		ct++;
	}
	else
	{
		if (IsFree)
		{
			// Evolve position
			x = x + dt*(vb+VXSPH) + 0.5*dt*dt*a;

			// Evolve velocity
			Vec3_t temp;
			temp = v;
			v = v + dt*a;
			vb = temp;

			// Evolve density
			if (ShepardFilter)
			{
				if (!isnan(SumDen/ZWab))
				{
					// Shepard filter
					Density = SumDen/ZWab;
					Densityb = Density;
				}
				else
				{
					Densityb = Density;
				}
			}
			else
			{
				double dens = Density;
				Density = Density + dt*dDensity;
				Densityb = dens;
			}

			// Evolve shear stress
			if (Material == 2)
			{
				Mat3_t RotationRateT;
				Mat3_t SRT,RS;
				Trans(Rotation,RotationRateT);
				Mult(ShearStress,RotationRateT,SRT);
				Mult(Rotation,ShearStress,RS);
				ShearStress = dt*(2.0*G*(StrainRate-1.0/(2.0+I(2,2))*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*I)+SRT+RS) + ShearStress;
				NormalStress = dt*(1.0/(2.0+I(2,2))*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*K) + NormalStress;
			}
		}
		ct=0;
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
