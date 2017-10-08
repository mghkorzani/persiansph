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

#ifndef SPH_DOMAIN_H
#define SPH_DOMAIN_H

#include <stdio.h>    // for NULL
#include <algorithm>  // for min,max

#include <hdf5.h>
#include <hdf5_hl.h>

#include <omp.h>

#include "Particle.h"
#include "Functions.h"
#include "Boundary_Condition.h"

namespace SPH {

class Domain
{
public:
	typedef void (*PtVel) (Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry);
	typedef void (*PtOut) (Particle * Particles, double & Prop1, double & Prop2,  double & Prop3);
	typedef void (*PtDom) (Domain & dom);
    // Constructor
    Domain();

    // Destructor
    ~Domain ();

    // Domain Part
    void AddSingleParticle			(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);		///< Add one particle
    void AddBoxLength				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz,double r, double Density,
							double h,int type, int rotation, bool random, bool Fixed);			///< Add a cube of particles with a defined dimensions
    void AddBoxNo				(int tag, Vec3_t const &V, size_t nx, size_t ny, size_t nz,double r, double Density,
							double h,int type, int rotation, bool random, bool Fixed);			///< Add a cube of particles with a defined numbers
    void DelParticles				(int const & Tags);									///< Delete particles by tag
    void CheckParticleLeave			();											///< Check if any particles leave the domain, they will be deleted

    void YZPlaneCellsNeighbourSearch		(int q1);										///< Create pairs of particles in cells of XZ plan
    void MainNeighbourSearch			();											///< Create pairs of particles in the whole domain
    void StartAcceleration			(Vec3_t const & a = Vec3_t(0.0,0.0,0.0));						///< Add a fixed acceleration such as the Gravity
    void PrimaryComputeAcceleration		();											///< Compute the solid boundary properties
    void LastComputeAcceleration		();											///< Compute the acceleration due to the other particles
    void CalcForce11				(Particle * P1, Particle * P2);								///< Calculates the contact force between fluid-fluid particles
    void CalcForce2233				(Particle * P1, Particle * P2);								///< Calculates the contact force between soil-soil/solid-solid particles
    void CalcForce12				(Particle * P1, Particle * P2);								///< Calculates the contact force between fluid-solid particles
    void CalcForce13				(Particle * P1, Particle * P2);								///< Calculates the contact force between fluid-soil particles
    void Move					(double dt);										///< Move particles

    void Solve					(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function

    void CellInitiate				();											///< Find the size of the domain as a cube, make cells and HOCs
    void ListGenerate				();											///< Generate linked-list
    void CellReset				();											///< Reset HOCs and particles' LL to initial value of -1

    void WriteXDMF				(char const * FileKey);									///< Save a XDMF file for the visualization
    void Save					(char const * FileKey);									///< Save the domain in a file
    void Load					(char const * FileKey);									///< Load the domain from the saved file
    void LoadResults				(char const * FileKey, double density);							///< Load the domain from one of the results

    void PrintInput				(char const * FileKey);
    void InFlowBCLeave				();
    void InFlowBCFresh				();
    void WholeVelocity				();
    void InitialChecks				();

    // Data
    Array <Particle*>				Particles; 	///< Array of particles
    double					R;		///< Particle Radius in addrandombox

    double					Time;    	///< The simulation time at each step
    double					AutoSaveInt;	///< Automatic save interval time step
    double					deltat;		///< Time Step
    double					deltatmin;	///< Minimum Time Step
    double					deltatint;	///< Initial Time Step
    double					DtAcc;		///< Coefficient for determining Time Step based on acceleration

    int 					Dimension;    	///< Dimension of the problem

    double					MuMax;		///< Max Dynamic viscosity for calculating the timestep
    double					CsMax;		///< Max speed of sound for calculating the timestep

    Vec3_t					Gravity;       	///< Gravity acceleration


    Vec3_t                 			TRPR;		///< Top right-hand point at rear of the domain as a cube
    Vec3_t                  			BLPF;           ///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                  			CellSize;      	///< Calculated cell size according to (cell size >= 2h)
    int		                		CellNo[3];      ///< No. of cells for linked list
    double 					Cellfac;	///< Factor which should be multiplied by h to change the size of cells (Min 2)
    double 					hmax;		///< Max of h for the cell size  determination
    Vec3_t                 			DomSize;	///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
    double					rhomax;

    int						*** HOC;	///< Array of "Head of Chain" for each cell

    size_t					VisEq;		///< Selecting variable to choose an equation for viscosity
    size_t					KernelType;	///< Selecting variable to choose a kernel
    size_t					SWIType;	///< Selecting variable to choose Soil-Water Interaction type
    bool					FSI;		///< Selecting variable to choose Fluid-Structure Interaction

    double 					XSPH;		///< Velocity correction factor
    double 					InitialDist;	///< Initial distance of particles for Inflow BC

    double					AvgVelocity;	///< Average velocity of the last two column for x periodic constant velocity

    size_t					Nproc;		///< No of threads which are going to use in parallel calculation
    omp_lock_t 					dom_lock;	///< Open MP lock to lock Interactions array
    Boundary					BC;
    PtOut					UserOutput;
    PtVel 					InCon;
    PtVel 					OutCon;
    PtVel 					AllCon;
    Vec3_t					DomMax;
    Vec3_t					DomMin;
    PtDom					GeneralBefore;	///< Pointer to a function: to modify particles properties before CalcForce function
    PtDom					GeneralAfter;	///< Pointer to a function: to modify particles properties after CalcForce function
    size_t					Scheme;		///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog
    bool					TimestepConstrain1;	///< Acceleration check for timestep size

    Array<Array<std::pair<size_t,size_t> > >	SMPairs;
    Array<Array<std::pair<size_t,size_t> > >	NSMPairs;
    Array<Array<std::pair<size_t,size_t> > >	FSMPairs;
    Array< size_t > 				FixedParticles;
    Array< size_t >				FreeFSIParticles;

    Array<std::pair<size_t,size_t> >		Initial;
    Mat3_t I;
    String					OutputName[3];
};

}; // namespace SPH

#include "Interaction.cpp"
#include "Domain.cpp"

#endif // SPH_DOMAIN_H
