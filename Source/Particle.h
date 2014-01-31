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

// MechSys
#include <mechsys/linalg/matvec.h>

namespace SPH {

class Particle
{
public:

    // Constructor
    Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double R0, double h0, bool Fixed=false);


    // Data
    bool   IsFree;                  ///< Check the particle if it is free to move or not
    Vec3_t  x;                       ///< Position of the particle
    Vec3_t  xb;                      ///< Previous position for Verlet algorithm
    Vec3_t  v;                       ///< Velocity of the particle
    Vec3_t  a;                       ///< Acceleration of the particle
    double Pressure;                ///< Pressure at the position of the particle
    double Density;                 ///< Density at the position of the particle
    double Densityb;                ///< Previous density for the Verlet algorithm
    double RefDensity;				 ///< Reference Density of Particle
    double Mass;                    ///< Mass of the particle
    double dDensity;                ///< Rate of density change in time
    double h;                       ///< Smoothing length of the particle
    double hr;                      ///< Reference smoothing length of the particle
    double R;                       ///< Radius of the particle
    int    ID;						 ///< an Integer value to identify type of the particles
    int    LL;						 ///< Linked-List variable to show next particle in list of a cell
    int    CC[3];					 ///< Current cell No for the particle (linked-list)


    // Methods
    void Move (double dt);                                                  ///< Update the important quantities of a particle
    bool CellUpdate  (Vec3_t CellSize, Vec3_t BLPF);                     	  ///< Check if the particle cell needs to be updated

};

inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double R0, double h0,bool Fixed)
{
    x = x0;
    xb = x;
    v = v0;
    Density = Density0;
    RefDensity = Density0;
    Densityb = Density;
    Mass = Mass0;
    IsFree = !Fixed;
    hr = h0;
    h = hr;
    R = R0;
    dDensity=0.0;
    Pressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    LL=0;
}

inline void Particle::Move (double dt)
{
    if (IsFree)
    {
        // Evolve position and velocity
        Vec3_t xa;
        xa = 2*x - xb + a*dt*dt;
        v = 0.5*(xa - xb)/dt;
        xb = x;
        x = xa;

        // Evolve density
        double dens = Density;
        Density = Densityb + 2*dt*dDensity;
        Densityb = dens;
    }
}

inline bool Particle::CellUpdate (Vec3_t CellSize, Vec3_t BLPF)
{
	bool update;
	//    if (CC == (int) (x(0) - BLPF(0)) / CellSize(0), (int) (x(1) - BLPF(1)) / CellSize(1), (int) (x(2) - BLPF(2)) / CellSize(2)) return false;
	//    else return true;
	    if ((CC[0] == (int) ((x(0) - BLPF(0)) / CellSize(0))) && (CC[1] == (int) ((x(1) - BLPF(1)) / CellSize(1))) && (CC[2] == 0)) update=false;
	    else update=true;
	    return update;
}

}; // namespace SPH

#endif // MECHSYS_SPH_PARTICLE_H
