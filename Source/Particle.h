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
    Vec3_t xo;                      ///< Previous position of the particle to check if it needs to move
    Vec3_t x;                       ///< Position of the particle
    Vec3_t xb;                      ///< Previous position for Verlet algorithm
    Vec3_t v;                       ///< Velocity of the particle
    Vec3_t a;                       ///< Acceleration of the particle
    double Pressure;                ///< Pressure at the position of the particle
    double Density;                 ///< Density at the position of the particle
    double Densityb;                ///< Previous density for the Verlet algorithm
    double RefDensity;				///< Reference Density of Particle
    double Mass;                    ///< Mass of the particle
    double dDensity;                ///< Rate of density change in time
    double h;                       ///< Smoothing length of the particle
    double hr;                      ///< Reference smoothing length of the particle
    double R;                       ///< Radius of the particle
    int    ID;						///< an Integer value to identify type of the particles


    // Methods
    void Move (double dt);                                                  ///< Update the important quantities of the simulation
    void StartAccel (Vec3_t acc = Vec3_t(0.0,0.0,0.0)) {a = acc;};          ///< Start the acceleration of the particle with one predefined value
    void Translate  (Vec3_t const & V) {x+=V; xb+=V;};                      ///< Translate a SPH particle a vector V
    void ResetDisplacements ();                                             ///< Reset the displacement for the Verlet algorithm
    double MaxDisplacement ();                                              ///< Find the maximum displacement

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
//  std::cout << Mass << std::endl;
    IsFree = !Fixed;
    hr = h0;
    h = hr;
    R = R0;
    dDensity=0.0;
    Pressure=0.0;
    ID = Tag;
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
    }
        // Evolve density
        double dens = Density;
        Density = Densityb + 2*dt*dDensity;
        Densityb = dens;


}

inline void Particle::ResetDisplacements ()
{
    xo = x;
}

inline double Particle::MaxDisplacement ()
{
    return Norm (xo-x);
}

}; // namespace SPH

#endif // MECHSYS_SPH_PARTICLE_H
