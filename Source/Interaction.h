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

#ifndef MECHSYS_SPH_INTERACTION_H
#define MECHSYS_SPH_INTERACTION_H


#include <Particle.h>
#include <Functions.h>

namespace SPH {

class Interaction
{
public:
    // Constructor
    Interaction (Particle * Pt1, Particle * Pt2); ///< Default constructor

    // Methods
    bool UpdateContacts (double alpha);          ///< Find neighbor particle and make contact
    void CalcForce      (double dt = 0.0);       ///< Calculates the contact force between particles


    // Data
    Particle * P1;                               ///< Pointer to first particle
    Particle * P2;                               ///< Pointer to second particle
    double alpha;                                ///< Coefficient of bulk viscosity
    double beta;                                 ///< Coefficient of Neumann-Richtmyer viscosity
    double h;                                    ///< Smoothing length
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Interaction::Interaction (Particle * Pt1, Particle * Pt2)
{
    P1 = Pt1;
    P2 = Pt2;
    h  = P1->h;                                  ///< It should be revised
    alpha = 0.75;
    beta = 2*alpha;
}

inline void Interaction::CalcForce(double dt)
{
    double di = P1->Density;
    double dj = P2->Density;
    double d0i = P1->Density0;
    double d0j = P2->Density0;
    Vec3_t vij = P2->v - P1->v;
    Vec3_t rij = P2->x - P1->x;
    double MUij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    double Cij = 0.5*(SoundSpeed(di)+SoundSpeed(dj));
    double PIij;
    if (dot(vij,rij)<0) PIij = (-alpha*Cij*MUij+beta*MUij*MUij)/(0.5*(di+dj));                           ///<(2.74) Li, Liu Book
    else                PIij = 0.0;
    P1->a += d0j*(Pressure(di)/(di*di)+Pressure(dj)/(dj*dj)+PIij)*rij*GradKernel(norm(rij),h)/norm(rij); ///<(2.73) Li, Liu Book
    P2->a -= d0i*(Pressure(di)/(di*di)+Pressure(dj)/(dj*dj)+PIij)*rij*GradKernel(norm(rij),h)/norm(rij);
    P1->dDensity += d0j*dot(vij,rij)*GradKernel(norm(rij),h)/norm(rij);                                  ///<(2.58) Li, Liu Book
    P2->dDensity += d0i*dot(vij,rij)*GradKernel(norm(rij),h)/norm(rij);
    P1->Pressure = Pressure(di);
    P2->Pressure = Pressure(dj);
}

inline bool Interaction::UpdateContacts (double alpha)
{
    if (Norm(P2->x-P1->x)<=P1->h+P2->h+2*alpha) return true;
    else return false;
}

}; // namespace SPH

#endif // MECHSYS_SPH_INTERACTION_H
