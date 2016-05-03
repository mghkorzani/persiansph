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

#ifndef MECHSYS_SPH_DOMAIN_H
#define MECHSYS_SPH_DOMAIN_H

// Std Lib
#include <stdio.h>			/// for NULL
#include <algorithm>		/// for min,max

#include "Particle.h"
#include "Functions.h"


// HDF File Output
#include <hdf5.h>
#include <hdf5_hl.h>


namespace SPH {

class Boundary
{
public:
    // Data
    double	inDensity;	///< Apply a certain density to inflow particles
    double	outDensity;	///< Apply a certain density to outflow particles
    double	allDensity;	///< Apply a certain density to outflow particles

    Vec3_t 	inv;		///< Apply a certain velocity to inflow particles
    Vec3_t 	outv;		///< Apply a certain velocity to outflow particles
    Vec3_t  allv;		///< Apply a certain velocity to all particle

    bool 	Periodic[3];///< Considering periodic in all directions => 0=X, 1=Y, 2=Z

    int 	InOutFlow;	///< Considering inflow in all directions  by adding and deleting particles=> [0]=X, [1]=Y, [2]=Z and 0=none, 1=-
    double  InFlowLoc1;
    double  InFlowLoc2;
    double  InFlowLoc3;
   double  OutFlowLoc;
    double	cellfac;
    int		inoutcounter;
    bool	MassConservation;

    Array <int>				OutPart;
	Array <int>				InPart;

	Boundary();
};

inline Boundary::Boundary()
{
	allv	= 0.0;
	inv		= 0.0;
	outv	= 0.0;
	Periodic[0]=Periodic[1]=Periodic[2] = false;
	inDensity = 0.0;
	outDensity = 0.0;
	allDensity = 0.0;
	InOutFlow = 0;
	InFlowLoc1 = 0.0;
	InFlowLoc2 = 0.0;
	InFlowLoc3 = 0.0;
	OutFlowLoc = 0.0;
	cellfac = 4.0;
	inoutcounter = 0;
	MassConservation = false;
}



class Domain
{
public:
	typedef void (*PtVel) (Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry);
	typedef void (*PtDom) (Domain & dom);
    // Constructor
    Domain();

    // Destructor
    ~Domain ();

    // Domain Part
    void AddSingleParticle			(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);						///< Add one particle
    void AddBoxLength				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz,double r, double Density, double h,
    									int type, int rotation, bool random, bool Fixed);												///< Add a cube of particles with a defined dimensions
    void AddBoxNo					(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double r, double Density, double h,
    									int type, int rotation, bool random, bool Fixed);												///< Add a cube of particles with a defined numbers

    void DelParticles				(int const & Tags);																					///< Delete particles by tag
    void CheckParticleLeave			();																									///< Check if any particles leave the domain, they will be deleted

    void YZPlaneCellsNeighbourSearch(int q1);																							///< Create pairs of particles in cells of XZ plan
    void MainNeighbourSearch		();																									///< Create pairs of particles in the whole domain
    void StartAcceleration			(Vec3_t const & a = Vec3_t(0.0,0.0,0.0));															///< Add a fixed acceleration such as the Gravity
    void PrimaryComputeAcceleration	();																									///< Compute the solid boundary properties
    void LastComputeAcceleration	();																									///< Compute the acceleration due to the other particles
    void CalcForce11				(Particle * P1, Particle * P2);																		///< Calculates the contact force between particles
    void CalcForce2233				(Particle * P1, Particle * P2);																		///< Calculates the contact force between particles
    void CalcForce13				(Particle * P1, Particle * P2);																		///< Calculates the contact force between particles
    void Move						(double dt);																						///< Move particles

    void Solve						(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);						///< The solving function

    void CellInitiate				();																									///< Find the size of the domain as a cube, make cells and HOCs
    void ListGenerate				();																									///< Generate linked-list
    void CellReset					();																									///< Reset HOCs and particles' LL to initial value of -1

    void WriteXDMF					(char const * FileKey);																				///< Save a XDMF file for the visualization
    void Save						(char const * FileKey);																				///< Save the domain in a file
    void Load						(char const * FileKey);																				///< Load the domain from the saved file
    void LoadResults				(char const * FileKey, double density);																///< Load the domain from one of the results

    void PrintInput					(char const * FileKey);
    void InFlowBCLeave				();
    void InFlowBCFresh				();
    void WholeVelocity				();
    void InitialChecks				();

    // Data
    Array <Particle*>		Particles;     	 	///< Array of particles
    double					R;					///< Particle Radius in addrandombox

    double					Time;          	 	///< The simulation time at each step
    double					AutoSaveInt;		///< Automatic save interval time step
    double					deltat;				///< Time Step

    int 					Dimension;    	  	///< Dimension of the problem

    double					MuMax;				///< Max Dynamic viscosity for calculating the timestep
    double					CsMax;				///< Max speed of sound for calculating the timestep

    Vec3_t					Gravity;       	 	///< Gravity acceleration


    Vec3_t                  TRPR;				///< Top right-hand point at rear of the domain as a cube
    Vec3_t                  BLPF;           	///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                  CellSize;      		///< Calculated cell size according to (cell size >= 2h)
    int		                CellNo[3];      	///< No. of cells for linked list
    double 					Cellfac;			///< Factor which should be multiplied by h to change the size of cells (Min 2)
	double 					hmax;				///< Max of h for the cell size  determination
    Vec3_t                  DomSize;			///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
	double					rhomax;

    int					*** HOC;				///< Array of "Head of Chain" for each cell

    size_t					VisEq;				///< Selecting variable to choose an equation for viscosity
    size_t					KernelType;			///< Selecting variable to choose a kernel
    size_t					SeepageType;		///< Selecting variable to choose a Seepage method

    double 					XSPH;				///< Velocity correction factor
    double 					InitialDist;		///< Initial distance of particles for calculation of tensile instability

    double					AvgVelocity;		///< Average velocity of the last two column for x periodic constant velocity

    size_t					Nproc;				///< No of threads which are going to use in parallel calculation
    omp_lock_t 				dom_lock;			///< Open MP lock to lock Interactions array
    Boundary				BC;
    PtVel 					InCon;
    PtVel 					OutCon;
    PtVel 					AllCon;
    Vec3_t					DomMax;
    Vec3_t					DomMin;
    PtDom					GeneralBefore;		///< Pointer to a function: to modify particles properties before CalcForce function
    PtDom					GeneralAfter;		///< Pointer to a function: to modify particles properties after CalcForce function
    size_t					Scheme;				///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog
    bool					TimestepConstrain1; ///< Acceleration check for timestep size

    Array<Array<std::pair<size_t,size_t> > >	Pairs;
    Array< size_t > 							FixedParticles;
    Array<std::pair<size_t,size_t> >			FixedPairs;
    Mat3_t I;
    int						Excemption;
};
/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////
void General(Domain & dom)
{
}

void InFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.inv;
	Den = bdry.inDensity;
}

void OutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.outv;
	Den = bdry.outDensity;
}

void AllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.allv;
	Den = bdry.allDensity;
}

// Constructor
inline Domain::Domain ()
{
    Excemption = 0;
    Time    = 0.0;
    AutoSaveInt = 0.0;

    Dimension = 2;
    DomSize	= 0.0,0.0,0.0;

    Gravity	= 0.0,0.0,0.0;
    SeepageType = 0;

    Cellfac = 2.0;

    KernelType	= 0;
    VisEq	= 0;
    Scheme	= 0;


    XSPH	= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;

    TRPR = 0.0;
    BLPF = 0.0;

    InCon = & InFlowCon;
    OutCon = & OutFlowCon;
    AllCon = & AllFlowCon;
    GeneralBefore = & General;
    GeneralAfter = & General;

    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
    I = OrthoSys::I;
    TimestepConstrain1 = true;
}

inline Domain::~Domain ()
{
	size_t Max = Particles.Size();
	for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

inline void Domain::AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
{
   	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,h,Fixed));
}

inline void Domain::AddBoxNo(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
{
    if (!(type==0 || type==1))
    {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90))
    {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 =>" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxNo with defined numbers of particles--------------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y;

    double qin = 0.03;
    srand(100);

    if (Dimension==3)
    {
   		double z;

    	if (type==0)
    	{
    		//Hexagonal close packing
 		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
				y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2.0*i+1)*r;
				y = V(1) + (2.0*j+1)*r;
				z = V(2) + (2.0*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		#pragma omp parallel for num_threads(Nproc)
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			Particles[i]->Mass = Mass;
		}
    }

    if (Dimension==2)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
    		    for (size_t j=0; j<ny; j++)
    		    for (size_t i=0; i<nx; i++)
    		    {
					x = V(0) + (2*i+(j%2)+1)*r;
					y = V(1) + (sqrt(3.0)*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
			}
    		else
    		{
    		    for (size_t i=0; i<nx; i++)
    		    for (size_t j=0; j<ny; j++)
    		    {
					x = V(0) + (sqrt(3.0)*i+1)*r;
					y = V(1) + (2*j+(i%2)+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2*i+1)*r;
				y = V(1) + (2*j+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
			}

    	}
    }

	R = r;
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
{
    if (!(type==0 || type==1))
    {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90))
    {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 => Cubic Close Packing" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y,xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);

    if (Dimension==3)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r))
			{
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
						y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
						z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						i++;
						if ((k%2!=0) && (j%2!=0)) xp = V(0) + (2*i+(j%2)+(k%2)-1)*r; else xp = V(0) + (2*i+(j%2)+(k%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				}
				k++;
				zp = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
			}
    	}
    	else
    	{
    		//Cubic packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r))
			{
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2.0*i+1)*r;
						y = V(1) + (2.0*j+1)*r;
						z = V(2) + (2.0*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+1)*r;
					}
					j++;
					yp = V(1) + (2.0*j+1)*r;
				}
				k++;
				zp = V(2) + (2.0*k+1)*r;
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		#pragma omp parallel for num_threads(Nproc)
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			Particles[i]->Mass = Mass;
		}
    }

    if (Dimension==2)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
				j = 0;
				yp = V(1);

				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2*i+(j%2)+1)*r;
						y = V(1) + (sqrt(3.0)*j+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+(j%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*j+1)*r;
				}
			}
    		else
    		{
				i = 0;
				xp = V(0);

				while (xp <= (V(0)+Lx-r))
				{
					j = 0;
					yp = V(1);
					while (yp <= (V(1)+Ly-r))
					{
						x = V(0) + (sqrt(3.0)*i+1)*r;
						y = V(1) + (2*j+(i%2)+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						j++;
						yp = V(1) + (2*j+(i%2)+1)*r;
					}
					i++;
					xp = V(0) + (sqrt(3.0)*i+1)*r;
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
    		j = 0;
			yp = V(1);

			while (yp <= (V(1)+Ly-r))
			{
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2*i+1)*r;
					y = V(1) + (2*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
					i++;
					xp = V(0) + (2*i+1)*r;
				}
				j++;
				yp = V(1) + (2*j+1)*r;
			}

    	}
    }

	R = r;
}

inline void Domain::DelParticles (int const & Tags)
{
    Array<int> idxs; // indices to be deleted

	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->ID==Tags)
		{
			omp_set_lock(&dom_lock);
        	idxs.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particles to delete");
    Particles.DelItems (idxs);

    std::cout << "\n" << "Particle(s) with Tag No. " << Tags << " has been deleted" << std::endl;
}

inline void Domain::CheckParticleLeave ()
{
	Array <int> DelParticles;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
    {
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			omp_set_lock(&dom_lock);
			DelParticles.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }

	if (DelParticles.Size()>0)
	{
		std::cout<< DelParticles.Size()<< " particle(s) left the Domain"<<std::endl;
		Particles.DelItems(DelParticles);
	}
}

inline void Domain::CellInitiate ()
{
	if (!(norm(TRPR)>0.0) && !(norm(BLPF)>0.0))
	{
		// Calculate Domain Size
		BLPF = Particles[0]->x;
		TRPR = Particles[0]->x;
		hmax = Particles[0]->h;
		rhomax = Particles[0]->Density;

		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > TRPR(0)) TRPR(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > TRPR(1)) TRPR(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > TRPR(2)) TRPR(2) = Particles[i]->x(2);

			if (Particles[i]->x(0) < BLPF(0)) BLPF(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) < BLPF(1)) BLPF(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) < BLPF(2)) BLPF(2) = Particles[i]->x(2);

			if (Particles[i]->h > hmax) hmax=Particles[i]->h;
			if (Particles[i]->Density > rhomax) rhomax=Particles[i]->Density;
			if (Particles[i]->Mu > MuMax) MuMax=Particles[i]->Mu;
			if (Particles[i]->Cs > CsMax) CsMax=Particles[i]->Cs;
		}
	}

	// Override the calculated domain size
	if (DomMax(0)>TRPR(0)) TRPR(0) = DomMax(0);
	if (DomMax(1)>TRPR(1)) TRPR(1) = DomMax(1);
	if (DomMax(2)>TRPR(2)) TRPR(2) = DomMax(2);
	if (DomMin(0)<BLPF(0)) BLPF(0) = DomMin(0);
	if (DomMin(1)<BLPF(1)) BLPF(1) = DomMin(1);
	if (DomMin(2)<BLPF(2)) BLPF(2) = DomMin(2);


	//Because of Hexagonal close packing in x direction domain is modified
	if (!BC.Periodic[0]) {TRPR(0) += hmax/2;	BLPF(0) -= hmax/2;}else{TRPR(0) += R; BLPF(0) -= R;}
	if (!BC.Periodic[1]) {TRPR(1) += hmax/2;	BLPF(1) -= hmax/2;}else{TRPR(1) += R; BLPF(1) -= R;}
	if (!BC.Periodic[2]) {TRPR(2) += hmax/2;	BLPF(2) -= hmax/2;}else{TRPR(2) += R; BLPF(2) -= R;}

    // Calculate Cells Properties
	switch (Dimension)
	{case 2:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		CellNo[2] = 1;

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],0.0);
		break;

	case 3:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(2)-BLPF(2))/(Cellfac*hmax)))-((TRPR(2)-BLPF(2))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[2] = int(ceil((TRPR(2)-BLPF(2))/(Cellfac*hmax)));
		else
			CellNo[2] = int(floor((TRPR(2)-BLPF(2))/(Cellfac*hmax)));

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	// Periodic BC modifications
	if (BC.Periodic[0]) CellNo[0] += 2;
    if (BC.Periodic[1]) CellNo[1] += 2;
    if (BC.Periodic[2]) CellNo[2] += 2;

    if (BC.Periodic[0]) DomSize[0] = (TRPR(0)-BLPF(0));
    if (BC.Periodic[1]) DomSize[1] = (TRPR(1)-BLPF(1));
    if (BC.Periodic[2]) DomSize[2] = (TRPR(2)-BLPF(2));

    // Initiate Head of Chain array for Linked-List
    HOC = new int**[(int) CellNo[0]];
    for(int i =0; i<CellNo[0]; i++){
       HOC[i] = new int*[CellNo[1]];
       for(int j =0; j<CellNo[1]; j++){
           HOC[i][j] = new int[CellNo[2]];
           for(int k = 0; k<CellNo[2];k++){
              HOC[i][j][k] = -1;
           }
       }
    }
}

inline void Domain::ListGenerate ()
{
	int i, j, k, temp=0;
	switch (Dimension)
	{case 2:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));

			if (i<0)
            {
                    if ((BLPF(0) - Particles[a]->x(0)) <= hmax) i=0;
                            else std::cout<<"Leaving i<0"<<std::endl;
            }
            if (j<0)
            {
                    if ((BLPF(1) - Particles[a]->x(1)) <= hmax) j=0;
                            else std::cout<<"Leaving j<0"<<std::endl;
            }
			if (i>=CellNo[0])
			{
					if ((Particles[a]->x(0) - TRPR(0)) <= hmax) i=CellNo[0]-1;
							else std::cout<<"Leaving i>=CellNo"<<std::endl;
			}
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1)) <= hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving j>=CellNo"<<std::endl;
            }

			temp = HOC[i][j][0];
			HOC[i][j][0] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = 0;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	case 3:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));
			k= (int) (floor((Particles[a]->x(2) - BLPF(2)) / CellSize(2)));

            if (i<0)
            {
                    if ((BLPF(0) - Particles[a]->x(0))<=hmax) i=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (j<0)
            {
                    if ((BLPF(1) - Particles[a]->x(1))<=hmax) j=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (k<0)
            {
                    if ((BLPF(2) - Particles[a]->x(2))<=hmax) k=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
			if (i>=CellNo[0])
			{
					if ((Particles[a]->x(0) - TRPR(0))<=hmax) i=CellNo[0]-1;
							else std::cout<<"Leaving"<<std::endl;
			}
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1))<=hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (k>=CellNo[2])
            {
                    if ((Particles[a]->x(2) - TRPR(2))<=hmax) k=CellNo[2]-1;
                            else std::cout<<"Leaving"<<std::endl;
            }

            temp = HOC[i][j][k];
			HOC[i][j][k] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = k;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	if (BC.Periodic[0])
	{
	   for(int j =0; j<CellNo[1]; j++)
	   for(int k =0; k<CellNo[2]; k++)
	   {
		  HOC[CellNo[0]-1][j][k] =  HOC[1][j][k];
		  HOC[CellNo[0]-2][j][k] =  HOC[0][j][k];
	   }
	}
	if (BC.Periodic[1])
	{
	   for(int i =0; i<CellNo[0]; i++)
	   for(int k =0; k<CellNo[2]; k++)
	   {
		  HOC[i][CellNo[1]-1][k] =  HOC[i][1][k];
		  HOC[i][CellNo[1]-2][k] =  HOC[i][0][k];
	   }
	}
	if (BC.Periodic[2])
	{
	   for(int i =0; i<CellNo[0]; i++)
	   for(int j =0; j<CellNo[1]; j++)
	   {
		  HOC[i][j][CellNo[2]-1] =  HOC[i][j][1];
		  HOC[i][j][CellNo[2]-2] =  HOC[i][j][0];
	   }
	}
}

inline void Domain::CellReset ()
{

	#pragma omp parallel for schedule (static) num_threads(Nproc)
    for(int i =0; i<CellNo[0]; i++)
    {
		for(int j =0; j<CellNo[1]; j++)
		for(int k =0; k<CellNo[2];k++)
		{
			HOC[i][j][k] = -1;
		}
    }

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t a=0; a<Particles.Size(); a++)
    {
    	Particles[a]->LL = -1;
    }

    FixedParticles.Clear();
}

inline void Domain::MainNeighbourSearch()
{
    int q1;

    if (BC.Periodic[0])
    {
		#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
		for (q1=1;q1<(CellNo[0]-1); q1++)	YZPlaneCellsNeighbourSearch(q1);
	}
    else
    {
		#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
    	for (q1=0;q1<CellNo[0]; q1++)	YZPlaneCellsNeighbourSearch(q1);
    }
}

inline void Domain::YZPlaneCellsNeighbourSearch(int q1)
{
	int q3,q2;
	Array<std::pair<size_t,size_t> > LocalPairs;
	Array<std::pair<size_t,size_t> > FixedLocalPairs;

	for (BC.Periodic[2] ? q3=1 : q3=0;BC.Periodic[2] ? (q3<(CellNo[2]-1)) : (q3<CellNo[2]); q3++)
	for (BC.Periodic[1] ? q2=1 : q2=0;BC.Periodic[1] ? (q2<(CellNo[1]-1)) : (q2<CellNo[1]); q2++)
	{
		if (HOC[q1][q2][q3]==-1) continue;
		else
		{
			int temp1, temp2;
			temp1 = HOC[q1][q2][q3];

			while (temp1 != -1)
			{
				// The current cell  => self cell interactions
				temp2 = Particles[temp1]->LL;
				while (temp2 != -1)
				{
					if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
						LocalPairs.Push(std::make_pair(temp1, temp2));
					temp2 = Particles[temp2]->LL;
				}

				// (q1 + 1, q2 , q3)
				if (q1+1< CellNo[0])
				{
					temp2 = HOC[q1+1][q2][q3];
					while (temp2 != -1)
					{
						if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
							LocalPairs.Push(std::make_pair(temp1, temp2));
						temp2 = Particles[temp2]->LL;
					}
				}

				// (q1 + a, q2 + 1, q3) & a[-1,1]
				if (q2+1< CellNo[1])
				{
					for (int i = q1-1; i <= q1+1; i++)
					{
						if (i<CellNo[0] && i>=0)
						{
							temp2 = HOC[i][q2+1][q3];
							while (temp2 != -1)
							{
								if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
									LocalPairs.Push(std::make_pair(temp1, temp2));
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}

				// (q1 + a, q2 + b, q3 + 1) & a,b[-1,1] => all 9 cells above the current cell
				if (q3+1< CellNo[2])
				{
					for (int j=q2-1; j<=q2+1; j++)
					for (int i=q1-1; i<=q1+1; i++)
					{
						if (i<CellNo[0] && i>=0 && j<CellNo[1] && j>=0)
						{
							temp2 = HOC[i][j][q3+1];
							while (temp2 != -1)
							{
								if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
									LocalPairs.Push(std::make_pair(temp1, temp2));
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}
				temp1 = Particles[temp1]->LL;
			}
		}
	}

	//Transferring LocalPairs array in each thread to Pairs (a global array)
	omp_set_lock(&dom_lock);
	Pairs.Push(LocalPairs);
	omp_unset_lock(&dom_lock);

}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
    {
    	if (Particles[i]->IsFree)
    	{
    		if (Particles[i]->Material < 3 && Particles[i]->FirstStep)
    			Particles[i]->Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);

        	// Tensile Instability for all soil and solid particles
            if (Particles[i]->Material > 1 && Particles[i]->TI > 0.0)
            {
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2)
				{
					double teta, Sigmaxx, Sigmayy, C, S;

					if ((Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1))!=0.0) teta = 0.5*atan(2.0*Particles[i]->Sigma(0,1)/(Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1))); else teta = M_PI/4.0;
					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[i]->Sigma(0,0) + 2.0*C*S*Particles[i]->Sigma(0,1) + S*S*Particles[i]->Sigma(1,1);
					Sigmayy = S*S*Particles[i]->Sigma(0,0) - 2.0*C*S*Particles[i]->Sigma(0,1) + C*C*Particles[i]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[i]->TI * Sigmaxx/(Particles[i]->Density*Particles[i]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[i]->TI * Sigmayy/(Particles[i]->Density*Particles[i]->Density); else Sigmayy = 0.0;
					Particles[i]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[i]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[i]->TIR(0,1) = Particles[i]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else
				{
					Mat3_t Vec,Val,VecT,temp;

					Rotation(Particles[i]->Sigma,Vec,VecT,Val);
					if (Val(0,0)>0) Val(0,0) = -Particles[i]->TI * Val(0,0)/(Particles[i]->Density*Particles[i]->Density); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -Particles[i]->TI * Val(1,1)/(Particles[i]->Density*Particles[i]->Density); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -Particles[i]->TI * Val(2,2)/(Particles[i]->Density*Particles[i]->Density); else Val(2,2) = 0.0;
					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[i]->TIR);
				}
            }
    	}
    	else
    	{
       		// Reset the pressure and the induced velocity for solid boundaries
    		Particles[i]->vb = 0.0;
    		Particles[i]->Pressure = 0.0;
            set_to_zero(Particles[i]->Sigma);
            set_to_zero(Particles[i]->ShearStress);
    	}



    	//Reset to zero for all particles
    	Particles[i]->a			= a;
    	Particles[i]->SatCheck	= false;
        Particles[i]->dDensity	= 0.0;
        Particles[i]->VXSPH		= 0.0;
        Particles[i]->ZWab		= 0.0;
        Particles[i]->SumDen	= 0.0;
        Particles[i]->Vis		= 0.0;
        Particles[i]->SumKernel	= 0.0;
        set_to_zero(Particles[i]->StrainRate);
        set_to_zero(Particles[i]->RotationRate);
    }
}

inline void Domain::PrimaryComputeAcceleration ()
{
	if (FixedParticles.Size()>0)
	{
		#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
		for (size_t k=0; k<Pairs.Size();k++)
		{
			for (size_t a=0; a<Pairs[k].Size();a++)
			{
				if (Particles[Pairs[k][a].first]->Material == Particles[Pairs[k][a].second]->Material)
				{
					if (!Particles[Pairs[k][a].first]->IsFree)
					{
						size_t i	= Pairs[k][a].first;
						size_t j	= Pairs[k][a].second;
						Vec3_t xij	= Particles[i]->x-Particles[j]->x;
						double h	= (Particles[i]->h+Particles[j]->h)/2.0;

						// Correction of xij for Periodic BC
						if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(Particles[i]->CC[0]>Particles[j]->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
						if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(Particles[i]->CC[1]>Particles[j]->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
						if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(Particles[i]->CC[2]>Particles[j]->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

						double K = Kernel(Dimension, KernelType, norm(xij), h);

						omp_set_lock(&Particles[i]->my_lock);
							Particles[i]->SumKernel									+= K;
							if (Particles[i]->Material < 3)	Particles[i]->Pressure	+= Particles[j]->Pressure * K + dot(Gravity,xij)*Particles[j]->Density*K;
							if (Particles[i]->Material > 1)	Particles[i]->Sigma		=  Particles[i]->Sigma + K * Particles[j]->Sigma;
							Particles[i]->vb += Particles[j]->v * K;
						omp_unset_lock(&Particles[i]->my_lock);

	//					omp_set_lock(&dom_lock);
	//		        	FixedPairs.Push(Pairs[k][a]);
	//					omp_unset_lock(&dom_lock);
					}
					if (!Particles[Pairs[k][a].second]->IsFree)
					{
						size_t i	= Pairs[k][a].first;
						size_t j	= Pairs[k][a].second;
						Vec3_t xij	= Particles[i]->x-Particles[j]->x;
						double h	= (Particles[i]->h+Particles[j]->h)/2.0;

						// Correction of xij for Periodic BC
						if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(Particles[i]->CC[0]>Particles[j]->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
						if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(Particles[i]->CC[1]>Particles[j]->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
						if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(Particles[i]->CC[2]>Particles[j]->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

						double K = Kernel(Dimension, KernelType, norm(xij), h);
						omp_set_lock(&Particles[j]->my_lock);
							Particles[j]->SumKernel									+= K;
							if (Particles[j]->Material < 3)	Particles[j]->Pressure	+= Particles[i]->Pressure * K + dot(Gravity,xij)*Particles[i]->Density*K;
							if (Particles[j]->Material > 1)	Particles[j]->Sigma		=  Particles[j]->Sigma + K * Particles[i]->Sigma;
							Particles[j]->vb += Particles[i]->v * K;
						omp_unset_lock(&Particles[j]->my_lock);

	//					omp_set_lock(&dom_lock);
	//		        	FixedPairs.Push(Pairs[k][a]);
	//					omp_unset_lock(&dom_lock);
					}
				}
				else
				{
					if (Particles[Pairs[k][a].first]->Material == 3)
					{
						if (!Particles[Pairs[k][a].first]->SatCheck)
							if (Particles[Pairs[k][a].second]->CC[1] >= Particles[Pairs[k][a].first]->CC[1])
								if (Particles[Pairs[k][a].second]->x(1) >= Particles[Pairs[k][a].first]->x(1))
									Particles[Pairs[k][a].first]->SatCheck = true;
					}
					if (Particles[Pairs[k][a].second]->Material == 3)
					{
						if (!Particles[Pairs[k][a].second]->SatCheck)
							if (Particles[Pairs[k][a].first]->CC[1] >= Particles[Pairs[k][a].second]->CC[1])
								if (Particles[Pairs[k][a].first]->x(1) >= Particles[Pairs[k][a].second]->x(1))
									Particles[Pairs[k][a].second]->SatCheck = true;
					}
				}
			}
		}

		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->Material == 3)
			{
				if (Particles[i]->SatCheck && !Particles[i]->IsSat)
				{
					Particles[i]->Mass			= Particles[i]->V*(Particles[i]->RefDensity - Particles[i]->RhoF);
					Particles[i]->Density		= Particles[i]->Density - Particles[i]->RhoF;
					Particles[i]->Densityb		= Particles[i]->Densityb - Particles[i]->RhoF;
					Particles[i]->RefDensity	= Particles[i]->RefDensity - Particles[i]->RhoF;
					Particles[i]->IsSat			= true;
				}
				if (!Particles[i]->SatCheck && Particles[i]->IsSat)
				{
					Particles[i]->Mass			= Particles[i]->V*(Particles[i]->RefDensity + Particles[i]->RhoF);
					Particles[i]->Density		= Particles[i]->Density + Particles[i]->RhoF;
					Particles[i]->Densityb		= Particles[i]->Densityb + Particles[i]->RhoF;
					Particles[i]->RefDensity	= Particles[i]->RefDensity + Particles[i]->RhoF;
					Particles[i]->IsSat			= false;
				}
			}
		}


		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<FixedParticles.Size(); i++)
			if (Particles[FixedParticles[i]]->SumKernel!= 0.0)
			{
				size_t a = FixedParticles[i];
				if (Particles[a]->Material < 3)	Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
				if (Particles[a]->Material > 1) Particles[a]->Sigma		= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
				Particles[a]->vb		= Particles[a]->vb/Particles[a]->SumKernel;

				// Tensile Instability for fixed soil and solid particles
				if (Particles[a]->Material > 1 && Particles[a]->TI > 0.0)
				{
					// XY plane must be used, It is very slow in 3D
					if (Dimension == 2)
					{
						double teta, Sigmaxx, Sigmayy, C, S;

						if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0) teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))); else teta = M_PI/4.0;
						C = cos(teta);
						S = sin(teta);
						Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
						Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
						if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
						if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
						Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
						Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
						Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
					}
					else
					{
						Mat3_t Vec,Val,VecT,temp;

						Rotation(Particles[a]->Sigma,Vec,VecT,Val);
						if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
						if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
						if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
						Mult(Vec,Val,temp);
						Mult(temp,VecT,Particles[a]->TIR);
					}
				}
			}
	}
}

inline void Domain::LastComputeAcceleration ()
{
    #pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (size_t k=0; k<Pairs.Size();k++)
	{
		for (size_t i=0; i<Pairs[k].Size();i++)
		{
			switch (Particles[Pairs[k][i].first]->Material*Particles[Pairs[k][i].second]->Material)
			{case 1:
				CalcForce11(Particles[Pairs[k][i].first],Particles[Pairs[k][i].second]);
				break;
			case 4:
				CalcForce2233(Particles[Pairs[k][i].first],Particles[Pairs[k][i].second]);
				break;
			case 9:
				CalcForce2233(Particles[Pairs[k][i].first],Particles[Pairs[k][i].second]);
				break;
			case 3:
				CalcForce13(Particles[Pairs[k][i].first],Particles[Pairs[k][i].second]);
				break;
		   default:
				std::cout<<"Out of Interaction types"<<std::endl;
				abort();
				break;
			}
		}
	}
	Pairs.Clear();
//	FixedPairs.Clear();

//	//Min time step check based on the acceleration
	if (TimestepConstrain1)
	{
		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<Particles.Size(); i++)
			if (deltat > (0.25*sqrt(Particles[i]->h/norm(Particles[i]->a))) )
			{
				std::cout << "Please decrease the time step to"<< (0.25*sqrt(Particles[i]->h/norm(Particles[i]->a))) << std::endl;
				abort();
			}
	}
}

inline void Domain::Move (double dt)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if (Particles[i]->IsFree)
		{
			if (!Particles[i]->InOut)
				Particles[i]->Move(dt,DomSize,TRPR,BLPF,Scheme,I);
			else
				Particles[i]->x(0) += dt*Particles[i]->v(0);
		}

}

inline void Domain::InFlowBCLeave()
{
	size_t a;
	Array <int> DelPart,TempPart;
	Array<std::pair<Vec3_t,size_t> > AddPart;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			Particles[i]->InOut	= false;
			omp_set_lock(&dom_lock);
			DelPart.Push(i);
			omp_unset_lock(&dom_lock);
		}

	if (BC.InOutFlow==1 || BC.InOutFlow==3)
	{
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			if(Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc1)
			{
				Vec3_t temp1	= Particles[BC.InPart[i]]->x;
				Particles[BC.InPart[i]]->InOut	= false;
				temp1(0) -=  (BC.InFlowLoc3-BC.InFlowLoc2+InitialDist);
				AddPart.Push(std::make_pair(temp1,BC.InPart[i]));
				TempPart.Push(i);
			}
		BC.InPart.DelItems(TempPart);
		TempPart.Clear();
	}

	if (AddPart.Size() >= DelPart.Size())
	{
		for (size_t i=0 ; i<DelPart.Size() ; i++)
		{
			a = DelPart[i];
			Particles[a]->PresEq	= Particles[AddPart[i].second]->PresEq;
			Particles[a]->Alpha		= Particles[AddPart[i].second]->Alpha;
			Particles[a]->Beta		= Particles[AddPart[i].second]->Beta;
			Particles[a]->Mu		= Particles[AddPart[i].second]->Mu;
			Particles[a]->MuRef		= Particles[AddPart[i].second]->MuRef;
			Particles[a]->m			= Particles[AddPart[i].second]->m;
			Particles[a]->T0		= Particles[AddPart[i].second]->T0;
			Particles[a]->RefDensity= Particles[AddPart[i].second]->RefDensity;
			Particles[a]->Material	= 1;
			Particles[a]->Cs		= Particles[AddPart[i].second]->Cs;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->ID 		= Particles[AddPart[i].second]->ID;
			Particles[a]->Mass 		= Particles[AddPart[i].second]->Mass;
			Particles[a]->h			= Particles[AddPart[i].second]->h;
			Particles[a]->InOut		= true;
			BC.InPart.Push(a);
		}
		for (size_t i=DelPart.Size() ; i<AddPart.Size() ; i++)
		{
			Particles.Push(new Particle(Particles[AddPart[i].second]->ID,AddPart[i].first,Particles[AddPart[i].second]->v,Particles[AddPart[i].second]->Mass,Particles[AddPart[i].second]->RefDensity,Particles[AddPart[i].second]->h,false));
			a = Particles.Size()-1;
			Particles[a]->PresEq	= Particles[AddPart[i].second]->PresEq;
			Particles[a]->Alpha		= Particles[AddPart[i].second]->Alpha;
			Particles[a]->Beta		= Particles[AddPart[i].second]->Beta;
			Particles[a]->Mu		= Particles[AddPart[i].second]->Mu;
			Particles[a]->MuRef		= Particles[AddPart[i].second]->MuRef;
			Particles[a]->m			= Particles[AddPart[i].second]->m;
			Particles[a]->T0		= Particles[AddPart[i].second]->T0;
			Particles[a]->RefDensity= Particles[AddPart[i].second]->RefDensity;
			Particles[a]->Material	= 1;
			Particles[a]->Cs		= Particles[AddPart[i].second]->Cs;
			Particles[a]->InOut		= true;
			BC.InPart.Push(a);
		}
		DelPart.Clear();
		AddPart.Clear();
	}

	// this condition satisfies BC.InOutFlow==2 and ParticleLeaveDomain automatically.
	if (AddPart.Size() < DelPart.Size())
	{
		for (size_t i=0 ; i<AddPart.Size() ; i++)
		{
			a = DelPart[i];
			Particles[a]->PresEq	= Particles[AddPart[i].second]->PresEq;
			Particles[a]->Alpha		= Particles[AddPart[i].second]->Alpha;
			Particles[a]->Beta		= Particles[AddPart[i].second]->Beta;
			Particles[a]->Mu		= Particles[AddPart[i].second]->Mu;
			Particles[a]->MuRef		= Particles[AddPart[i].second]->MuRef;
			Particles[a]->m			= Particles[AddPart[i].second]->m;
			Particles[a]->T0		= Particles[AddPart[i].second]->T0;
			Particles[a]->RefDensity= Particles[AddPart[i].second]->RefDensity;
			Particles[a]->Material	= 1;
			Particles[a]->Cs		= Particles[AddPart[i].second]->Cs;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->ID 		= Particles[AddPart[i].second]->ID;
			Particles[a]->Mass 		= Particles[AddPart[i].second]->Mass;
			Particles[a]->h			= Particles[AddPart[i].second]->h;
			Particles[a]->InOut		= true;
//			BC.InPart.Push(a);
		}
		for (size_t i=AddPart.Size() ; i<DelPart.Size() ; i++)
		{
			TempPart.Push(DelPart[i]);
		}
		Particles.DelItems(TempPart);
		BC.inoutcounter = 1;
		DelPart.Clear();
		AddPart.Clear();
	}
}

inline void Domain::InFlowBCFresh()
{
	int temp, temp1;
	int q1,q2,q3;
	if (BC.inoutcounter == 0)
	{
		if (BC.InOutFlow==1 || BC.InOutFlow==3)
		{
			BC.InPart.Clear();
			BC.InFlowLoc1  = BLPF(0) + BC.cellfac*hmax;
			temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

			for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
			for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
			for (q1=0; q1<(temp1 + 1)                                      ; q1++)
			{
				if (HOC[q1][q2][q3]!=-1)
				{
					temp = HOC[q1][q2][q3];
					while (temp != -1)
					{
						if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) )
						{
							BC.InPart.Push(temp);
							Particles[temp]->InOut = true;
						}
						temp = Particles[temp]->LL;
					}
				}
			}
			BC.InFlowLoc2  = Particles[BC.InPart[0]]->x(0);
			BC.InFlowLoc3  = Particles[BC.InPart[0]]->x(0);
			#pragma omp parallel for schedule(static) num_threads(Nproc)
			for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			{
				if (Particles[BC.InPart[i]]->x(0) < BC.InFlowLoc2) BC.InFlowLoc2  = Particles[BC.InPart[i]]->x(0);
				if (Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc3) BC.InFlowLoc3  = Particles[BC.InPart[i]]->x(0);
			}
		}

		if (BC.InOutFlow==2 || BC.InOutFlow==3)
			BC.OutFlowLoc = TRPR(0) - BC.cellfac*hmax;

		BC.inoutcounter = 2;
	}

	if (BC.inoutcounter == 1)
	{
		BC.InPart.Clear();
		temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

		for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=0; q1<(temp1 + 1)                                      ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) && Particles[temp]->InOut)
						BC.InPart.Push(temp);
					temp = Particles[temp]->LL;
				}
			}
		}
		BC.inoutcounter = 2;
	}


	if (BC.InOutFlow==2 || BC.InOutFlow==3)
	{
		BC.OutPart.Clear();
		temp1 = (int) (floor((BC.OutFlowLoc - BLPF(0)) / CellSize(0)));

		for (q2=0     ; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0     ; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=temp1 ; q1<CellNo[0]                                        ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) >= BC.OutFlowLoc) )
					{
						BC.OutPart.Push(temp);
						Particles[temp]->InOut = true;
					}
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	Vec3_t vel;
	double den;

	if (BC.InPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
		{
			InCon(Particles[BC.InPart[i]]->x,vel,den,BC);
			Particles[BC.InPart[i]]->Material	= 1;
			if (norm(BC.inv)>0.0)
			{
				Particles[BC.InPart[i]]->v  = vel;
				Particles[BC.InPart[i]]->vb = vel;
			}
			if (BC.inDensity>0.0)
			{
				Particles[BC.InPart[i]]->Density  = den;
				Particles[BC.InPart[i]]->Densityb = den;
//				Particles[BC.InPart[i]]->RefDensity = BC.inDensity;
    			Particles[BC.InPart[i]]->Pressure = EOS(Particles[BC.InPart[i]]->PresEq, Particles[BC.InPart[i]]->Cs, Particles[BC.InPart[i]]->P0,Particles[BC.InPart[i]]->Density, Particles[BC.InPart[i]]->RefDensity);
			}
		}

	double temp11;
	if (BC.MassConservation) temp11 = (BC.InPart.Size()*1.0)/(BC.OutPart.Size()*1.0); else temp11 = 1.0;
//	std::cout<<temp11<<std::endl;

	if (BC.OutPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.OutPart.Size() ; i++)
		{
			OutCon(Particles[BC.OutPart[i]]->x,vel,den,BC);
			if (norm(BC.outv)>0.0)
			{
				Particles[BC.OutPart[i]]->v  = temp11*vel;
				Particles[BC.OutPart[i]]->vb = temp11*vel;
			}
			if (BC.outDensity>0.0)
			{
				Particles[BC.OutPart[i]]->Density  = den;
				Particles[BC.OutPart[i]]->Densityb = den;
//				Particles[BC.OutPart[i]]->RefDensity = BC.outDensity;
    			Particles[BC.OutPart[i]]->Pressure = EOS(Particles[BC.OutPart[i]]->PresEq, Particles[BC.OutPart[i]]->Cs, Particles[BC.OutPart[i]]->P0,Particles[BC.OutPart[i]]->Density, Particles[BC.OutPart[i]]->RefDensity);
			}
		}
}

inline void Domain::WholeVelocity()
{
    //Apply a constant velocity to all particles in the initial time step
    if (norm(BC.allv)>0.0 || BC.allDensity>0.0)
    {
    	Vec3_t vel = 0.0;
    	double den = 0.0;

		#pragma omp parallel for schedule (static) private(vel,den) num_threads(Nproc)
    	for (size_t i=0 ; i<Particles.Size() ; i++)
    	{
			AllCon(Particles[i]->x,vel,den,BC);
    		if (Particles[i]->IsFree && norm(BC.allv)>0.0)
    		{
				Particles[i]->v  = vel;
    		}
    		if (Particles[i]->IsFree && BC.allDensity>0.0)
    		{
				Particles[i]->Density  = den;
				Particles[i]->Densityb = den;
//				Particles[i]->RefDensity = BC.allDensity;
    			Particles[i]->Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
    		}
    	}
    }
}

inline void Domain::InitialChecks()
{
    if (KernelType==4) Cellfac = 3.0; else Cellfac = 2.0;
    if (Dimension == 2) I(2,2) = 0;

    if (BC.InOutFlow>0 && BC.Periodic[0])
    {
    	std::cout << "Periodic BC in the X direction cannot be used with In/Out-Flow BC simultaneously" << std::endl;
		abort();
    }

    if (Dimension<=1 || Dimension>3)
    {
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
    }
    if (KernelType==1)
    	if (VisEq>1)
    	{
        	std::cout << "Quadratic kernel can not be used with the second order viscosity equations" << std::endl;
    		abort();
    	}
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx)
{
//    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

    size_t idx_out = 1;
    double tout = Time;
    deltat = dt;

    size_t save_out = 1;
    double sout = AutoSaveInt;

    InitialChecks();
    CellInitiate();
    ListGenerate();
    PrintInput(TheFileKey);
    WholeVelocity();
//    std::cout<< Excemption<<std::endl;

    while (Time<tf && idx_out<=maxidx)
    {
    	StartAcceleration(Gravity);

    	if (BC.InOutFlow>0) InFlowBCFresh();

    	MainNeighbourSearch();

    	PrimaryComputeAcceleration();

    	GeneralBefore(*this);

    	LastComputeAcceleration();

    	GeneralAfter(*this);

        // output
        if (Time>=tout)
        {
            if (TheFileKey!=NULL)
           	{
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                WriteXDMF    (fn.CStr());
                std::cout << "\n" << "Output No. " << idx_out << " at " << Time << " has been generated" << std::endl;
           	}
            idx_out++;
            tout += dtOut;
        }
    	Move(dt);

        // Auto Save
       if (AutoSaveInt>0)
       {
    	   if (Time>=sout)
    	   {
    		   if (TheFileKey!=NULL)
    		   {
				   String fn;
				   fn.Printf    ("Auto Save_%s_%04d", TheFileKey, save_out);
				   Save   		(fn.CStr());
				   std::cout << "\n" << "Auto Save No. " << save_out << " at " << Time << " has been generated" << std::endl;
    		   }
    		save_out++;
    		sout += AutoSaveInt;
    	   }
       }

       Time += dt;

       if (BC.InOutFlow>0) InFlowBCLeave(); else CheckParticleLeave ();

       CellReset();

       ListGenerate();
    }
    std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

inline void Domain::PrintInput(char const * FileKey)
{
    //Writing Inputs in a Log file
	String fn(FileKey);
    std::ostringstream oss;

    oss << "Dimension = "<< Dimension << "D\n";

    oss << "Kernel Type = ";
	switch (KernelType)
    {
		case 0:
			oss << "0 => Qubic Spline\n";
			break;
		case 1:
			oss << "1 => Quadratic\n";
			break;
		case 2:
			oss << "2 => Quintic\n";
			break;
		case 3:
			oss << "3 => Gaussian with compact support of q<2\n";
			break;
		case 4:
			oss << "4 => Quintic Spline\n";
			break;
    }

//    oss << "Viscosity Equation = ";
//    if (Alpha!=0.0 || Beta!=0.0) oss << "Artificial Viscosity by Monaghan\nAlpha = "<<Alpha<<"   Beta = "<<Beta <<"\n";
//    else
//	{
//    	switch (VisEq)
//		{
//			case 0:
//				oss << "0 => Morris et al 1997\n";
//				break;
//			case 1:
//				oss << "1 => Shao et al 2003\n";
//				break;
//			case 2:
//				oss << "2 => Real viscosity for incompressible fluids\n";
//				break;
//			case 3:
//				oss << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)\n";
//				break;
//		}
//	}

//    oss << "Equation of State = ";
//	switch (PresEq)
//    {
//		case 0:
//			oss << "0 => P0+(Cs*Cs)*(Density-Density0)\n";
//			break;
//		case 1:
//			oss << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)\n";
//			break;
//		case 2:
//			oss << "2 => (Cs*Cs)*Density\n";
//			break;
//    }
//	oss << "Cs = "<<Cs<<" m/s\n";
//	oss << "P0 = "<<P0<<" Pa\n";
	oss << "\n";

    oss << "Min domain Point = " << BLPF <<" m\n";
    oss << "Max domain Point = " << TRPR <<" m\n";
    oss << "Max of smoothing length, h = " << hmax << " m\n";
    oss << "Cell factor in Linked List (at least 2) = " << Cellfac << "\n";
    oss << "Cell Size = " << CellSize <<" m\n";
    oss << "No of Cells in X Direction = " << CellNo[0] <<"\n" ;
    oss << "No of Cells in Y Direction = " << CellNo[1] <<"\n" ;
    oss << "No of Cells in Z Direction = " << CellNo[2] <<"\n" ;
	oss << "\n";

	oss << " Total No of Particles = " << Particles.Size() << "\n" ;

    // Check the time step
    double t1,t2;
    t1 = 0.25*hmax/(CsMax);
    if (MuMax>0.0) t2 = 0.125*hmax*hmax*rhomax/MuMax; else t2 =1000000.0;

    oss << "Max time step should be less than Min value of { "<< t1 <<" , "<< t2 <<" } S\n";
    oss << "Time Step = "<<deltat << " S\n";

    if (deltat > std::min(t1,t2))
    {
        std::cout << "Max time step should be less than Min value of { "<< t1 <<" , "<< t2 <<" }" << std::endl;
        std::cout << "Time Step = "<<deltat << std::endl;
    	std::cout << "Please decrease the time step and run again"<< std::endl;
    	abort();
    }
    oss << "External Acceleration = "<<Gravity<< " m/s2\n";
    oss << "No of Thread = "<<Nproc<<"\n";
//    oss << "Shepard Filter for Density = " << (Shepard ? "True" : "False") << "\n";
    oss << "Periodic Boundary Condition X dir= " << (BC.Periodic[0] ? "True" : "False") << "\n";
    oss << "Periodic Boundary Condition Y dir= " << (BC.Periodic[1] ? "True" : "False") << "\n";
    oss << "Periodic Boundary Condition Z dir= " << (BC.Periodic[2] ? "True" : "False") << "\n";

    fn = FileKey;
    fn.append("log.dat");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::WriteXDMF (char const * FileKey)
{
	String fn(FileKey);
    fn.append(".hdf5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    float * Posvec   = new float[3*Particles.Size()];
    float * Velvec   = new float[3*Particles.Size()];
    float * ACCvec   = new float[3*Particles.Size()];
    float * Pressure = new float[  Particles.Size()];
    float * ShearRate= new float[  Particles.Size()];
    float * Density  = new float[  Particles.Size()];
    float * Mass	 = new float[  Particles.Size()];
    float * sh	     = new float[  Particles.Size()];
    int   * Tag      = new int  [  Particles.Size()];
    int   * IsFree   = new int  [  Particles.Size()];
    float * Sigma    = new float[6*Particles.Size()];
    float * Strain   = new float[6*Particles.Size()];


    for (size_t i=0;i<Particles.Size();i++)
    {
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        ACCvec  [3*i  ] = float(Particles[i]->a(0));
        ACCvec  [3*i+1] = float(Particles[i]->a(1));
        ACCvec  [3*i+2] = float(Particles[i]->a(2));
       	Pressure[i    ] = float(Particles[i]->Pressure);
        ShearRate[i   ] = float(Particles[i]->ShearRate);
        Density [i    ] = float(Particles[i]->Density);
        Mass	[i    ] = float(Particles[i]->Mass);
        sh	    [i    ] = float(Particles[i]->h);
        Tag     [i    ] = int  (Particles[i]->ID);
        if (Particles[i]->IsFree)
        	IsFree[i] = int  (1);
        else
        	IsFree[i] = int  (0);
        Sigma  [6*i  ] = float(Particles[i]->Sigma(0,0));
        Sigma  [6*i+1] = float(Particles[i]->Sigma(0,1));
        Sigma  [6*i+2] = float(Particles[i]->Sigma(0,2));
        Sigma  [6*i+3] = float(Particles[i]->Sigma(1,1));
        Sigma  [6*i+4] = float(Particles[i]->Sigma(1,2));
        Sigma  [6*i+5] = float(Particles[i]->Sigma(2,2));
        Strain [6*i  ] = float(Particles[i]->Strain(0,0));
        Strain [6*i+1] = float(Particles[i]->Strain(0,1));
        Strain [6*i+2] = float(Particles[i]->Strain(0,2));
        Strain [6*i+3] = float(Particles[i]->Strain(1,1));
        Strain [6*i+4] = float(Particles[i]->Strain(1,2));
        Strain [6*i+5] = float(Particles[i]->Strain(2,2));
   }

    int data[1];
    String dsname;
    hsize_t dims[1];
    dims[0]=1;
    data[0]=Particles.Size();
    dsname.Printf("/NP");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,data);
    dims[0] = 3*Particles.Size();
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("Acceleration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ACCvec);
    dims[0] = Particles.Size();
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    dsname.Printf("ShearRate");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ShearRate);
    dsname.Printf("Density");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density);
    dsname.Printf("Mass");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Mass);
    dsname.Printf("h");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,sh);
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);
    dsname.Printf("IsFree");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,IsFree);
    dims[0] = 6*Particles.Size();
    dsname.Printf("Sigma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);



    delete [] Posvec;
    delete [] Velvec;
    delete [] ACCvec;
    delete [] Pressure;
    delete [] ShearRate;
    delete [] Density;
    delete [] Mass;
    delete [] sh;
    delete [] Tag;
    delete [] IsFree;

   //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"SPHCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"10\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Acceleration \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Position\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ShearRate\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ShearRate \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::Save (char const * FileKey)
{
    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");

    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing the number of particles in the domain
    int data[1];
    data[0]=Particles.Size();
    hsize_t dims[1];
    dims[0]=1;
    H5LTmake_dataset_int(file_id,"/NP",1,dims,data);

    for (size_t i=0; i<Particles.Size(); i++)
    {
        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gcreate2(file_id, par.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Storing some scalar variables
        double dat[1];
        dat[0] = Particles[i]->Mass;
        H5LTmake_dataset_double(group_id,"Mass",1,dims,dat);
        dat[0] = Particles[i]->Density;
        H5LTmake_dataset_double(group_id,"Rho",1,dims,dat);
        dat[0] = Particles[i]->h;
        H5LTmake_dataset_double(group_id,"h",1,dims,dat);

        int tag[1];
        tag[0] = Particles[i]->ID;
        H5LTmake_dataset_int(group_id,"Tag",1,dims,tag);

        // Change  to integer for fixity
        int IsFree[1];
        if (Particles[i]->IsFree == true) IsFree[0]=1;
        		else IsFree[0]=0;
        H5LTmake_dataset_int(group_id,"IsFree",1,dims,IsFree);

        // Storing vectorial variables
        double cd[3];
        hsize_t dd[1];
        dd[0] = 3;

        cd[0]=Particles[i]->x(0);
        cd[1]=Particles[i]->x(1);
        cd[2]=Particles[i]->x(2);
        H5LTmake_dataset_double(group_id,"x",1,dd,cd);

        cd[0]=Particles[i]->v(0);
        cd[1]=Particles[i]->v(1);
        cd[2]=Particles[i]->v(2);
        H5LTmake_dataset_double(group_id,"v",1,dd,cd);
    }

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
}

inline void Domain::Load (char const * FileKey)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"/NP",data);
    size_t NP = data[0];

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gopen2(file_id, par.CStr(),H5P_DEFAULT);

        Particles.Push(new Particle(0,Vec3_t(0,0,0),Vec3_t(0,0,0),0,0,0,false));

        // Loading vectorial variables
        double cd[3];
        H5LTread_dataset_double(group_id,"x",cd);
        Particles[Particles.Size()-1]->x = Vec3_t(cd[0],cd[1],cd[2]);

//        H5LTread_dataset_double(group_id,"v",cd);
//        Particles[Particles.Size()-1]->v = Vec3_t(cd[0],cd[1],cd[2]);

        // Loading the scalar quantities of the particle
        double dat[1];
        H5LTread_dataset_double(group_id,"Mass",dat);
        Particles[Particles.Size()-1]->Mass = dat[0];

        H5LTread_dataset_double(group_id,"Rho",dat);
        Particles[Particles.Size()-1]->Density = dat[0];
        Particles[Particles.Size()-1]->RefDensity = dat[0];			// Because of the constructor in Particle
        Particles[Particles.Size()-1]->Densityb = dat[0];			// Because of the constructor in Particle

        H5LTread_dataset_double(group_id,"h",dat);
        Particles[Particles.Size()-1]->h = dat[0];			// Because of the constructor in Particle

        int datint[1];
        H5LTread_dataset_int(group_id,"Tag",datint);
        Particles[Particles.Size()-1]->ID = datint[0];

        H5LTread_dataset_int(group_id,"IsFree",datint);
        if (datint[0] == 1) Particles[Particles.Size()-1]->IsFree=true;
        		else Particles[Particles.Size()-1]->IsFree=false;
    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

inline void Domain::LoadResults (char const * FileKey, double density)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"NP",data);
    size_t NP = data[0];
//    std::cout<<NP<<std::endl;

    float * Posvec   = new float[3*data[0]];
    float * Velvec   = new float[3*data[0]];
    float * Pressure = new float[  data[0]];
    float * Density  = new float[  data[0]];
    float * Mass	 = new float[  data[0]];
    float * sh	     = new float[  data[0]];
    int   * Tag      = new int  [  data[0]];
    int   * IsFree   = new int  [  data[0]];

    H5LTread_dataset_float(file_id, "Position"	, Posvec);
    H5LTread_dataset_float(file_id, "Velocity"	, Velvec);
    H5LTread_dataset_float(file_id, "Pressure"	, Pressure);
    H5LTread_dataset_float(file_id, "Density"	, Density);
    H5LTread_dataset_float(file_id, "Mass"		, Mass);
    H5LTread_dataset_float(file_id, "h"			, sh);
    H5LTread_dataset_int  (file_id, "Tag"		, Tag);
    H5LTread_dataset_int  (file_id, "IsFree"	, IsFree);

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        Particles.Push(new Particle(Tag[i],Vec3_t(Posvec[3*i],Posvec[3*i+1],Posvec[3*i+2]),Vec3_t(Velvec[3*i],Velvec[3*i+1],Velvec[3*i+2]),Mass[i],Density[i],sh[i],!(IsFree[i])));

        Particles[Particles.Size()-1]->Pressure = Pressure[i];
        Particles[Particles.Size()-1]->RefDensity = density;			// Because of the constructor in Particle

    }

    delete [] Posvec;
    delete [] Velvec;
    delete [] Pressure;
    delete [] Density;
    delete [] Mass;
    delete [] sh;
    delete [] Tag;
    delete [] IsFree;

    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

}; // namespace SPH

#endif // MECHSYS_SPH_DOMAIN_H
