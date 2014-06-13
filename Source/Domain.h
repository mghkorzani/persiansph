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


#include <Interaction.h>

// HDF File Output
#include <hdf5.h>
#include <hdf5_hl.h>

#include <mechsys/util/stopwatch.h>
#include <mechsys/util/string.h>

namespace SPH {



class Domain
{
public:

    // Constructor
    Domain();

    // Destructor
    ~Domain ();

    // Methods
    void AddSingleParticle			(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);						///< Add one particle
    void AddBoxLength				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    									double Mass, double Density, double h, bool Fixed);												///< Add a cube of particles with a defined dimensions
    void AddRandomBox				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    									double Mass, double Density, double h, size_t RandomSeed=100);									///< Add a cube of random positioned particles with a defined dimensions

    void DelParticles				(int const & Tags);																					///< Delete particles by tag
    void CheckParticleLeave			();																									///< Check if any particle are leaving the domain or not

    void StartAcceleration			(Vec3_t const & a = Vec3_t(0.0,0.0,0.0));															///< Add a fixed acceleration such as the Gravity
    void ComputeAcceleration		(double dt);																						///< Compute the acceleration due to the other particles
    void Move						(double dt);																						///< Move particles

    void ConstVel					();																									///< Give a defined velocity to the first and last 4 column of cells in X direction
    void ConstVelPart2				();																									///< Give a zero acceleration to all constant velocity zone (Correction for the first part)
    void AvgParticleVelocity		();																									///< Calculate the average velocity of whole domain (free particles)
    void Solve						(double tf, double dt, double dtOut, char const * TheFileKey, size_t Nproc);						///< The solving function

    void InitiateInteractions		();																									///< Reset the interaction array and make the first Interaction and PInteraction array
    void CellInitiate				();																									///< Find the size of the domain as a cube, make cells and HOCs
    void ListGenerate				();																									///< Generate linked-list

    void NeighbourSearch			(int q1, int q2, int q3);																			///< Find neighbor cells and make interaction for cell(q1,q2,q3)
    void CreateInteraction			(int a, int b);																						///< Create interaction between two particles and put it in PInteractions

    void CellReset					();																									///< Reset HOCs to initial value of -1
    void ListandInteractionUpdate	();																									///< Checks if linked-list needs to be updated
    void UpdateInteractions			();																									///< Update PInteraction (possible interaction) list

    void WriteXDMF					(char const * FileKey);																				///< Save a XDMF file for the visualization
    void Save						(char const * FileKey);																				///< Save the domain in a file
    void Load						(char const * FileKey);																				///< Load the domain from a file

    // Data
    Array <Particle*>		Particles;     	 	///< Array of particles
    Array <Interaction*>	Interactions;  	 	///< Array of interactions
    Array <Interaction*>	PInteractions; 	 	///< Array of possible interactions

    Array <int>				IinSI;				///< Array to save Interaction No for fixed particles

    double					Time;          	 	///< The simulation time at each step
    double					AutoSaveInt;		///< Automatic save interval

    int 					Dimension;    	  	///< Dimension of the problem

    double 					Alpha;				///< Artificial Viscosity Alpha Factor
    double					Beta;				///< Artificial Viscosity Beta Factor
    double					MU;					///< Dynamic Viscosity

    Vec3_t					Gravity;       	 	///< Gravity acceleration
    double					Cs;					///< Spead of sound

    Vec3_t                  TRPR;				///< Top right-hand point at rear of the domain as a cube
    Vec3_t                  BLPF;           	///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                  CellSize;      		///< Calculated cell size according to (cell size >= 2h)
    int		                CellNo[3];      	///< No. of cells for linked list
    double 					Cellfac;			///< Factor which should be multiplied by h to change the size of cells (Min 2)
	double 					hmax;				///< Max of h for the cell size  determination

    int					*** HOC;				///< Array of "Head of Chain" for each cell
    int					 ** ExInteract;			///< Array to save existing interaction No. of "Interactions" to use in "PInteractions"

    bool					RigidBody; 			///< If it is true, A rigid body with RTag will be considered
    int						RBTag;				///< Tag of particles for rigid body
    Vec3_t					RBForce;			///< Rigid body force

    double 					P0;					///< background pressure for equation of state
    int						PresEq;				///< Selecting variable to choose a equation of state

    bool					Periodic; 			///< If it is true, periodic boundary condition along x direction will be considered
    double					ConstVelPeriodic;	///< Define a constant velocity in the left and the right side of the domain in x direction in periodic condition

    bool					PressureBoundary;	///< If it is true, it will get max pressure for solid boundary from neighbors but if false, mean value.

    double 					XSPH;				///< Velocity correction factor
    double 					TI;					///< Tensile instability factor
    double 					InitialDist;		///< Initial distance of particles for calculation of tensile instability

    double					AvgVelocity;		///< Average velocity of the whole domain

};
/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Constructor
inline Domain::Domain ()
{
    Time    = 0.0;
    AutoSaveInt = 0.0;

    Dimension = 2;

    Alpha	= 0.0;
    Beta	= 0.0;
    MU		= 0.0;

    Gravity	= 0.0,0.0,0.0;
    Cs		= 0.0;

    Cellfac = 2.0;

    RigidBody = false;
    RBTag	= 0;
    RBForce = 0.0;

    P0		= 0.0;
    PresEq	= 0;

    Periodic = false;
    ConstVelPeriodic = 0.0;

    PressureBoundary = false;

    XSPH	= 0.0;
    TI		= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;
}

inline Domain::~Domain ()
{
	size_t Max = Interactions.Size();
	for (size_t i=1; i<=Max ; i++) Interactions.DelItem(Max-i);

	Max = Particles.Size();
	for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

// Methods
inline void Domain::AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
{
   	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,h,Fixed));
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, bool Fixed)
{
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length--------------" << std::endl;

    size_t PrePS = Particles.Size();

    for (size_t i=0; i<nx; i++)
    for (size_t j=0; j<ny; j++)
    for (size_t k=0; k<nz; k++)
    {
		double x = V(0)+i*Lx/nx;
		double y = V(1)+j*Ly/ny;
		double z = V(2)+k*Lz/nz;
		Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,h,Fixed));
    }

    std::cout << "\nNo. of added particles = " << Particles.Size()-PrePS << std::endl;
}

inline void Domain::AddRandomBox(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, size_t RandomSeed)
{
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating random packing of particles by AddRandomBox----------------" << std::endl;

    size_t PrePS = Particles.Size();

    double qin = 0.95;
    srand(RandomSeed);
    for (size_t i=0; i<nx; i++)
    for (size_t j=0; j<ny; j++)
    for (size_t k=0; k<nz; k++)
    {
		double x = V(0)+(i+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lx/nx;
		double y = V(1)+(j+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Ly/ny;
		double z = V(2)+(k+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lz/nz;
		Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,h,false));
    }

    std::cout << "\nNo. of added particles = " << Particles.Size()-PrePS << std::endl;
}

inline void Domain::DelParticles (int const & Tags)
{
    Array<int> idxs; // indices to be deleted

	#pragma omp parallel for
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->ID==Tags) idxs.Push(i);
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particle to be deleted");
    Particles.DelItems (idxs);

    std::cout << "\n" << "Particles with Tag No. " << Tags << " has been deleted" << std::endl;
}

inline void Domain::CheckParticleLeave ()
{
	Array <int> DelParticles;

	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
    {
		if (!Periodic){
			if ((Particles[i]->x(0) >= TRPR(0)) || (Particles[i]->x(1) >= TRPR(1)) || (Particles[i]->x(2) >= TRPR(2)) ||
					(Particles[i]->x(0) <= BLPF(0)) || (Particles[i]->x(1) <= BLPF(1)) || (Particles[i]->x(2) <= BLPF(2))) DelParticles.Push(i);}
		else{
			if ((Particles[i]->x(1) >= TRPR(1)) || (Particles[i]->x(2) >= TRPR(2)) ||
					(Particles[i]->x(1) <= BLPF(1)) || (Particles[i]->x(2) <= BLPF(2))) DelParticles.Push(i);}
    }

	if (DelParticles.Size()>0)
	{
		std::cout<< DelParticles.Size()<< " Particles are leaving the Domain"<<std::endl;

		for(size_t i =0; i<Particles.Size(); i++)
		for(size_t j =0; j<Particles.Size(); j++)
		{
			ExInteract[i][j] = -1;
		}

		Particles.DelItems(DelParticles);

		size_t Max = Interactions.Size();
		for (size_t i=1; i<=Max ; i++) Interactions.DelItem(Max-i);
		PInteractions.Resize(0);
	    Interactions.Resize(0);

		CellReset();
	    for (size_t a=0; a<Particles.Size(); a++)
	    	{
	        Particles[a]->LL = -1;
	    	}
	    ListGenerate();
	    UpdateInteractions();
	}
}

inline void Domain::CellInitiate ()
{
	// Calculate Domain Size
	BLPF = Particles[0]->x;
	TRPR = Particles[0]->x;

	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
    {
		if (Particles[i]->x(0) > TRPR(0)) TRPR(0) = Particles[i]->x(0);
        if (Particles[i]->x(1) > TRPR(1)) TRPR(1) = Particles[i]->x(1);
        if (Particles[i]->x(2) > TRPR(2)) TRPR(2) = Particles[i]->x(2);

        if (Particles[i]->x(0) < BLPF(0)) BLPF(0) = Particles[i]->x(0);
        if (Particles[i]->x(1) < BLPF(1)) BLPF(1) = Particles[i]->x(1);
        if (Particles[i]->x(2) < BLPF(2)) BLPF(2) = Particles[i]->x(2);

        if (Particles[i]->h > hmax) hmax=Particles[i]->h;
    }

	TRPR += hmax;
	BLPF -= hmax;

    if ((BLPF(0) < 0.0) | (BLPF(0) < 0.0) | (BLPF(0) < 0.0))
    	{
    	std::cout << "Problem to allocate Cells!" << std::endl;
    	std::cout << "Particle with minus coordinate" << std::endl;
    	abort();
    	}

    // Calculate Cells Properties
	switch (Dimension)
	{case 1:
		CellNo[0] = (int) (floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		CellNo[1] = 1;
		CellNo[2] = 1;
		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],0.0,0.0);
		break;
	case 2:
		CellNo[0] = (int) (floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		CellNo[1] = (int) (floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		CellNo[2] = 1;
		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],0.0);
		break;
	case 3:
		CellNo[0] = (int) (floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		CellNo[1] = (int) (floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		CellNo[2] = (int) (floor((TRPR(2)-BLPF(2))/(Cellfac*hmax)));
		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);
		break;
	default:
		std::cout << "Please correct dimension and run again" << std::endl;
		break;
	}

    if (Periodic) CellNo[0] += 2;

    std::cout << "No of Cells in X Direction = " << CellNo[0] << std::endl;
    std::cout << "No of Cells in Y Direction = " << CellNo[1] << std::endl;
    std::cout << "No of Cells in Z Direction = " << CellNo[2] << std::endl;

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
	{case 1:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			temp = HOC[i][0][0];
			HOC[i][0][0] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = 0;
			Particles[a]->CC[2] = 0;
		}
		break;
	case 2:
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
            if (Periodic)
            {
            	if (i>=CellNo[0]-2)
				{
						if ((Particles[a]->x(0) - TRPR(0)) <= hmax) i=CellNo[0]-3;
								else std::cout<<"Leaving i>=CellNo"<<std::endl;
				}
            }
            else
            {
            	if (i>=CellNo[0])
				{
						if ((Particles[a]->x(0) - TRPR(0)) <= hmax) i=CellNo[0]-1;
								else std::cout<<"Leaving i>=CellNo"<<std::endl;
				}

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
            if (Periodic)
            {
            	if (i>=CellNo[0])
				{
						if ((Particles[a]->x(0) - TRPR(0))<=hmax) i=CellNo[0]-3;
								else std::cout<<"Leaving"<<std::endl;
				}
            }
            else
            {
            	if (i>=CellNo[0])
				{
						if ((Particles[a]->x(0) - TRPR(0))<=hmax) i=CellNo[0]-1;
								else std::cout<<"Leaving"<<std::endl;
				}

            }
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1))<=hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving"<<std::endl;
            }

            temp = HOC[i][j][k];
			HOC[i][j][k] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = k;
		}
		break;
	}

	if (Periodic)
	{
	   for(int j =0; j<CellNo[1]; j++)
	   for(int k =0; k<CellNo[2];k++)
	   {
		  HOC[CellNo[0]-1][j][k] =  HOC[1][j][k];
		  HOC[CellNo[0]-2][j][k] =  HOC[0][j][k];
	   }
	}
}

inline void Domain::CellReset ()
{
    for(int i =0; i<CellNo[0]; i++)
    for(int j =0; j<CellNo[1]; j++)
    for(int k =0; k<CellNo[2];k++)
    {
    	HOC[i][j][k] = -1;
    }
}

inline void Domain::ListandInteractionUpdate()
{
	for (size_t i=0; i<Particles.Size(); i++)
    {
	if (Particles[i]->CellUpdate(CellSize,BLPF))
		{
		CellReset();
	    for (size_t a=0; a<Particles.Size(); a++)
	    	{
	        Particles[a]->LL = -1;
	    	}
	    ListGenerate();
	    UpdateInteractions();
	    break;
		}
    }
}

inline void Domain::CreateInteraction(int a, int b)
{
	if (Particles[a]->IsFree || Particles[b]->IsFree)
    {
		if (ExInteract [a][b] == -1)
		{
			Interactions.Push(new Interaction(Particles[a],Particles[b],Dimension,Alpha,Beta,Cs,MU,XSPH,TI,InitialDist,P0,PresEq));
			ExInteract [a][b] = Interactions.size() - 1;
			ExInteract [b][a] = Interactions.size() - 1;
			PInteractions.Push(Interactions[ ExInteract [a][b] ]);
		}
		else
		{
			PInteractions.Push(Interactions[ ExInteract [a][b] ]);
		}

		// Pressure boundary
		if (!Particles[a]->IsFree || !Particles[b]->IsFree)
		{
			IinSI.Push(ExInteract [a][b]);
		}
	}
}

inline void Domain::NeighbourSearch(int q1, int q2, int q3)
{
	int temp1, temp2;
	temp1 = HOC[q1][q2][q3];

	while (temp1 != -1)
	{
		// The current cell  => selfcell interactions
		temp2 = Particles[temp1]->LL;
		while (temp2 != -1)
		{
			CreateInteraction(temp1, temp2);
			temp2 = Particles[temp2]->LL;
		}

		// (q1 + 1, q2 , q3)
		if (q1+1< CellNo[0])
		{
			temp2 = HOC[q1+1][q2][q3];
			while (temp2 != -1)
			{
				CreateInteraction(temp1, temp2);
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
						CreateInteraction(temp1, temp2);
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
						CreateInteraction(temp1, temp2);
						temp2 = Particles[temp2]->LL;
					}
				}
			}
		}
		temp1 = Particles[temp1]->LL;
	}
}

inline void Domain::InitiateInteractions()
{
	// Initiate ExInteract
	ExInteract = new int*[Particles.Size()];
    for(size_t i =0; i<Particles.Size(); i++)
    	{
    	ExInteract[i] = new int[Particles.Size()];
    	for(size_t j =0; j<Particles.Size(); j++)
    		{
    		ExInteract[i][j] = -1;
    		}
    	}

	// Delete old interactions
	size_t Max = Interactions.Size();
	for (size_t i=1; i<=Max ; i++) Interactions.DelItem(Max-i);
    PInteractions.Resize(0);
    IinSI.Resize(0);

    if (!Periodic)
    {
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=0; q1<CellNo[0]; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else  NeighbourSearch(q1,q2,q3);
		}
    }
    else
    {
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=1; q1<CellNo[0]-1; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else  NeighbourSearch(q1,q2,q3);
		}
    }
}

inline void Domain::UpdateInteractions()
{
	PInteractions.Resize(0);
    IinSI.Resize(0);

    if (!Periodic)
    {
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=0; q1<CellNo[0]; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else  NeighbourSearch(q1,q2,q3);
		}
    }
    else
    {
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=1; q1<CellNo[0]-1; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else  NeighbourSearch(q1,q2,q3);
		}
    }
}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->a = a;
        Particles[i]->dDensity = 0.0;
        Particles[i]->VXSPH=0.0;
    }
}

inline void Domain::ComputeAcceleration (double dt)
{

	if (PressureBoundary)
	{
		// Get max value of pressure from neighbor for solid body particles
		#pragma omp parallel for
		for (size_t i=0; i<Particles.Size(); i++) if (!Particles[i]->IsFree) Particles[i]->Pressure=0.0;

		#pragma omp parallel for
		for (size_t i=0; i<IinSI.Size(); i++)
		{
			if (Interactions[IinSI [i]]->P1->IsFree && !Interactions[IinSI [i]]->P2->IsFree)
			{
				if (Interactions[IinSI [i]]->P2->Pressure < Interactions[IinSI [i]]->P1->Pressure) Interactions[IinSI [i]]->P2->Pressure = Interactions[IinSI [i]]->P1->Pressure;
			}
			else
			{
				if (Interactions[IinSI [i]]->P1->Pressure < Interactions[IinSI [i]]->P2->Pressure) Interactions[IinSI [i]]->P1->Pressure = Interactions[IinSI [i]]->P2->Pressure;

			}
		}
	}
	else
	{
		// Get max value of pressure from neighbor for solid body particles
		#pragma omp parallel for
		for (size_t i=0; i<Particles.Size(); i++) if (!Particles[i]->IsFree) Particles[i]->Pressure=0.0;

		#pragma omp parallel for
		for (size_t i=0; i<IinSI.Size(); i++)
		{
			if (Interactions[IinSI [i]]->P1->IsFree && !Interactions[IinSI [i]]->P2->IsFree)
			{
				if (Interactions[IinSI [i]]->P2->Pressure == 0.0) Interactions[IinSI [i]]->P2->Pressure = Interactions[IinSI [i]]->P1->Pressure;
				else Interactions[IinSI [i]]->P2->Pressure = (Interactions[IinSI [i]]->P2->Pressure + Interactions[IinSI [i]]->P1->Pressure)/2;
			}
			else
			{
				if (Interactions[IinSI [i]]->P1->Pressure == 0.0) Interactions[IinSI [i]]->P1->Pressure = Interactions[IinSI [i]]->P2->Pressure;
				else Interactions[IinSI [i]]->P1->Pressure = (Interactions[IinSI [i]]->P1->Pressure + Interactions[IinSI [i]]->P2->Pressure)/2;
			}
		}
	}

	#pragma omp parallel for
	for (size_t i=0; i<PInteractions.Size(); i++) PInteractions[i]->CalcForce(dt);
}

inline void Domain::Move (double dt)
{
	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt,Periodic,TRPR(0),BLPF(0),hmax);

	if (RigidBody)
	{
		RBForce = 0.0;
	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
		{
		if (Particles[i]->ID==RBTag) RBForce +=Particles[i]->Mass*Particles[i]->a;
		}
	}
}

inline void Domain::ConstVel ()
{
	if (Periodic && ConstVelPeriodic>0.0)
	{
		int temp;
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=CellNo[0]-6; q1<CellNo[0]-2; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree) Particles[temp]->v = Vec3_t (ConstVelPeriodic,0.0,0.0);
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	if (Periodic && ConstVelPeriodic>0.0)
	{
		int temp;
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=0; q1<4; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree) Particles[temp]->v = Vec3_t (ConstVelPeriodic,0.0,0.0);
					temp = Particles[temp]->LL;
				}
			}
		}
	}
}

inline void Domain::ConstVelPart2 ()
{
	if (Periodic && ConstVelPeriodic>0.0)
	{
		int temp;
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=CellNo[0]-6; q1<CellNo[0]-2; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree) Particles[temp]->a = Vec3_t (0.0,0.0,0.0);
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	if (Periodic && ConstVelPeriodic>0.0)
	{
		int temp;
		for (int q3=0; q3<CellNo[2]; q3++)
		for (int q2=0; q2<CellNo[1]; q2++)
		for (int q1=0; q1<4; q1++)
		{
			if (HOC[q1][q2][q3]==-1) continue;
			else
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree) Particles[temp]->a = Vec3_t (0.0,0.0,0.0);
					temp = Particles[temp]->LL;
				}
			}
		}
	}
}

inline void Domain::AvgParticleVelocity ()
{
	AvgVelocity = 0.0;
	int j = 0;

	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
		{
		if (Particles[i]->IsFree)
			{
			AvgVelocity += Particles[i]->v(0);
			j++;
			}
		}

	AvgVelocity = AvgVelocity / j;
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t Nproc)
{
    std::cout << "\nTotal No. of particles = " << Particles.Size()<< std::endl;
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

    size_t idx_out = 1;
    double tout = Time;

    size_t save_out = 1;
    double sout = AutoSaveInt;

    CellInitiate();
    ListGenerate();
    InitiateInteractions();

    // Initialization of constant velocity
	if (Periodic && ConstVelPeriodic>0.0)
	{
		#pragma omp parallel for
		for (size_t i=0; i<Particles.Size(); i++) if (Particles[i]->IsFree) Particles[i]->v = Vec3_t (ConstVelPeriodic,0.0,0.0);
	}

    while (Time<tf)
    {
    	StartAcceleration(Gravity);
    	ConstVel();
    	ComputeAcceleration(dt);
    	ConstVelPart2();
    	Move(dt);
    	AvgParticleVelocity ();

        // output
        if (Time>=tout)
        {
            if (TheFileKey!=NULL)
            	{
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                WriteXDMF    (fn.CStr());
                std::cout << "\n" << "Output No. " << idx_out << " at " << Time << " has been generated" << std::endl;
                std::cout << AvgVelocity<< std::endl;
            	}
            idx_out++;
            tout += dtOut;
        }

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

       CheckParticleLeave();
       ListandInteractionUpdate();
    }
}

inline void Domain::WriteXDMF (char const * FileKey)
{
	String fn(FileKey);
    fn.append(".hdf5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    float * Posvec   = new float[3*Particles.Size()];
    float * Velvec   = new float[3*Particles.Size()];
    float * Pressure = new float[  Particles.Size()];
    float * Density  = new float[  Particles.Size()];
    float * Mass	 = new float[  Particles.Size()];
    float * sh	     = new float[  Particles.Size()];
    int   * Tag      = new int  [  Particles.Size()];
    int   * IsFree   = new int  [  Particles.Size()];
    float * Force    = new float[3];


    for (size_t i=0;i<Particles.Size();i++)
    {
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        Pressure[i    ] = float(Particles[i]->Pressure);
        Density [i    ] = float(Particles[i]->Density);
        Mass	[i    ] = float(Particles[i]->Mass);
        sh	    [i    ] = float(Particles[i]->h);
        Tag     [i    ] = int  (Particles[i]->ID);
        if (Particles[i]->IsFree)
        	IsFree[i] = int  (1);
        else
        	IsFree[i] = int  (0);
    }
    Force  [0] = float(RBForce(0));
    Force  [1] = float(RBForce(1));
    Force  [2] = float(RBForce(2));

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
    dims[0] = Particles.Size();
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
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
    dims[0] = 3;
    dsname.Printf("Rigid_Body_Force");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Force);



    delete [] Posvec;
    delete [] Velvec;
    delete [] Pressure;
    delete [] Density;
    delete [] Mass;
    delete [] sh;
    delete [] Tag;
    delete [] IsFree;
    delete [] Force;

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
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
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
        dat[0] = Particles[i]->hr;
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
        Particles[Particles.Size()-1]->hr = dat[0];
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

}; // namespace SPH

#endif // MECHSYS_SPH_DOMAIN_H
