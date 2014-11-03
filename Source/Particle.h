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
#include <mechsys/linalg/matvec.h>

namespace SPH {

class Particle
{
public:

    // Constructor
    Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed=false);


    // Data
    bool   	IsFree;			///< Check the particle if it is free to move or not
    int    	ID;				///< an Integer value to identify type of the particles

    Vec3_t  x;				///< Position of the particle n
    Vec3_t  vb;				///< Velocity of the particle n-1
    Vec3_t  v;				///< Velocity of the particle n+1,
    Vec3_t	VXSPH;			///< Mean Velocity of neighbor particles
    Vec3_t  a;				///< Acceleration of the particle
    Vec3_t  Vis;			///< To check viscosity force

    double	ZWab;			///< Summation of the mb/db*Wab
    double	SumDen;			///< Summation of neighbor particle density
    double 	Pressure;		///< Pressure at the position of the particle
    double	Density;		///< Density at the position of the particle n+1
    double 	Densityb;		///< Density at the position of the particle n-1
    double 	RefDensity;		///< Reference Density of Particle
    double 	dDensity;		///< Rate of density change in time
    double 	Mass;			///< Mass of the particle

    double 	h;				///< Smoothing length of the particle
    double 	hr;				///< Reference smoothing length of the particle

    int    	LL;				///< Linked-List variable to show the next particle in the list of a cell
    int    	CC[3];			///< Current cell No for the particle (linked-list)
    Vec3_t 	NoSlip1;		///< Current cell No for the particle (linked-list)
    Vec3_t 	NoSlip2;		///< Current cell No for the particle (linked-list)

    int		ct;				///< Correction step for the Verlet Algorithm and Shepard filter

    omp_lock_t my_lock;		///< Open MP lock

    // Methods
    void Move			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, bool ShepardFilter);	///< Update the important quantities of a particle
    void translate		(Vec3_t acceleration, Vec3_t vellim, double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
};

inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0,bool Fixed)
{
	ct =0;
    x = x0;
    vb = v = v0;
    Densityb = Density = Density0;
    RefDensity = Density0;
    Mass = Mass0;
    IsFree = !Fixed;
    hr = h0;
    h = hr;
    dDensity=0.0;
    Pressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    NoSlip1=0.0;
    NoSlip2=0.0,0.0,1000000.0;
    LL=0;
    omp_init_lock(&my_lock);
    VXSPH = 0.0;
    ZWab = 0.0;
}

inline void Particle::Move (double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, bool ShepardFilter)
{
	if (ct<30)
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
		}
		else
		{
			// Evolve density for boundary particles
			double dens = Density;
			Density = Densityb + 2*dt*dDensity;
			Densityb = dens;
		}

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
			if (ShepardFilter && !isnan(SumDen/ZWab))
			{
				// Shepard filter
				Density = SumDen/ZWab;
				Densityb = Density;
			}
			else
			{
				double dens = Density;
				Density = Density + dt*dDensity;
				Densityb = dens;
			}
		}
		else
		{
			// Evolve density for boundary particles
			if (ShepardFilter && !isnan(SumDen/ZWab))
			{
				// Shepard filter
				Density = SumDen/ZWab;
				Densityb = Density;
			}
			else
			{
				double dens = Density;
				Density = Density + dt*dDensity;
				Densityb = dens;
			}
		}

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

		ct=0;
	}
}

inline void Particle::translate	(Vec3_t acceleration, Vec3_t vellim, double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
{
	if (v(0) < vellim(0) || v(1) < vellim(1) || v(2) < vellim(2))
	{
		x = x + dt*vb + 0.5*dt*dt*acceleration;

		Vec3_t temp;
		temp = v;
		v = vb + 2*dt*acceleration;
		vb = temp;
	}
	else
	{
		x = x + dt*v;
	}

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
