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
//#include <Force.hpp>


// HDF File Output
#include <hdf5.h>
#include <hdf5_hl.h>

//#include <mechsys/util/stopwatch.h>
//#include <mechsys/util/string.h>

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
    double  OutFlowLoc;
    double	cellfac;
    int		inoutcounter;
    bool	MassConservation;

    Array <int>				OutPart;
	Array <int>				InPart;

	Boundary();
};

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
    void CalcForce					(Particle * P1, Particle * P2);																		///< Calculates the contact force between particles
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
    void InFlowBCReset				();
    void WholeVelocity				();
    void InitialChecks				();

    // Data
    Array <Particle*>		Particles;     	 	///< Array of particles
    double					R;					///< Particle Radius in addrandombox

    double					Time;          	 	///< The simulation time at each step
    double					AutoSaveInt;		///< Automatic save interval time step
    double					deltat;				///< Time Step

    int 					Dimension;    	  	///< Dimension of the problem

    double 					Alpha;				///< Artificial Viscosity Alpha Factor
    double					Beta;				///< Artificial Viscosity Beta Factor
    double					MuMax;				///< Max Dynamic viscosity for calculating the timestep

    Vec3_t					Gravity;       	 	///< Gravity acceleration
    double					Cs;					///< Speed of sound
    double 					P0;					///< background pressure for equation of state


    Vec3_t                  TRPR;				///< Top right-hand point at rear of the domain as a cube
    Vec3_t                  BLPF;           	///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                  CellSize;      		///< Calculated cell size according to (cell size >= 2h)
    int		                CellNo[3];      	///< No. of cells for linked list
    double 					Cellfac;			///< Factor which should be multiplied by h to change the size of cells (Min 2)
	double 					hmax;				///< Max of h for the cell size  determination
    Vec3_t                  DomSize;			///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
	double					rhomax;

    int					*** HOC;				///< Array of "Head of Chain" for each cell

    size_t					PresEq;				///< Selecting variable to choose an equation of state
    size_t					VisEq;				///< Selecting variable to choose an equation for viscosity
    size_t					KernelType;			///< Selecting variable to choose a kernel

    bool					NoSlip;				///< To simulate No Slip Boundary
    double 					XSPH;				///< Velocity correction factor
    double 					TI;					///< Tensile instability factor
    double 					InitialDist;		///< Initial distance of particles for calculation of tensile instability
    bool					Shepard;			///< It is a first order correction for the density which is called Shepard Filter

    double					AvgVelocity;		///< Average velocity of the last two column for x periodic constant velocity

    size_t					Nproc;				///< No of threads which are going to use in parallel calculation
    omp_lock_t 				dom_lock;			///< Open MP lock to lock Interactions array
    Boundary				BC;
    PtVel 					InCon;
    PtVel 					OutCon;
    PtVel 					AllCon;
    Vec3_t					DomMax;
    Vec3_t					DomMin;
    bool 					BCDensityUpdate;
    PtDom					GeneralBefore;		///< Pointer to a function: to modify particles properties before CalcForce function
    PtDom					GeneralAfter;		///< Pointer to a function: to modify particles properties after CalcForce function

    Array<Array<std::pair<size_t,size_t> > > Pairs;
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

inline Boundary::Boundary()
{
	allv	= 0.0;
	inv		= 0.0;
	outv	= 0.0;
	Periodic[0]=Periodic[1]=Periodic[2]			= false;
	inDensity = 0.0;
	outDensity = 0.0;
	allDensity = 0.0;
	InOutFlow = 0;
	InFlowLoc1 = 0.0;
	OutFlowLoc = 0.0;
	cellfac = 4.0;
	inoutcounter = 0;
	MassConservation = false;
}

// Constructor
inline Domain::Domain ()
{
    Time    = 0.0;
    AutoSaveInt = 0.0;

    Dimension = 2;
    DomSize	= 0.0,0.0,0.0;
    Alpha	= 0.0;
    Beta	= 0.0;

    Gravity	= 0.0,0.0,0.0;
    Cs		= 0.0;

    Cellfac = 2.0;


    P0		= 0.0;
    PresEq	= 0;
    KernelType	= 0;
    VisEq	= 0;

    NoSlip	= false;
    BCDensityUpdate = true;

    XSPH	= 0.0;
    TI		= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;
    Shepard = true;

    TRPR = 0.0;
    BLPF = 0.0;

    InCon = & InFlowCon;
    OutCon = & OutFlowCon;
    AllCon = & AllFlowCon;
    GeneralBefore = & General;
    GeneralAfter = & General;

    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
}

inline Domain::~Domain ()
{
	size_t Max = Particles.Size();
	for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

inline void Domain::CalcForce(Particle * P1, Particle * P2)
{
	double h = (P1->h+P2->h)/2;
	double di,dj;
    double mi = P1->Mass;
    double mj = P2->Mass;
    double Mui = P1->Mu;
    double Muj = P2->Mu;

    if (!P1->IsFree) di = DensitySolid(PresEq, Cs, P0,P1->Pressure, P2->RefDensity); else di = P1->Density;
    if (!P2->IsFree) dj = DensitySolid(PresEq, Cs, P0,P2->Pressure, P1->RefDensity); else dj = P2->Density;

    Vec3_t vij = P1->v - P2->v;
    Vec3_t xij = P1->x - P2->x;

    // Correction of xij for Periodic BC
    if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

    double rij	= norm(xij);
    double GK	= GradKernel(Dimension, KernelType, rij, h) / rij;
    double K	= Kernel(Dimension, KernelType, rij, h);

    //Artificial Viscosity
    double PIij = 0.0;
    if (Alpha!=0.0 || Beta!=0.0)
    {
    	double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    	double Cij = 0.5*(SoundSpeed(PresEq, Cs, di, P1->RefDensity)+SoundSpeed(PresEq, Cs, dj, P2->RefDensity));
    	if (dot(vij,xij)<0) PIij = (-Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
    	else                PIij = 0.0;
    }

    //Tensile Instability
    double TIij = 0.0, TIji = 0.0;
    if ((TI > 0.0) && (P1->Pressure < 0.0) && (P2->Pressure < 0.0))
    {
        TIij = TIji= TI*(-P1->Pressure/(di*di)-P2->Pressure/(dj*dj))*pow((K/Kernel(Dimension, KernelType, InitialDist, h)),4);
        if (!P1->IsFree) TIij = 0.0;
        if (!P2->IsFree) TIji = 0.0;
    }

	//Real Viscosity
    Mat3_t StrainRate;
    Vec3_t VI = 0.0;
    if (!NoSlip || (P1->IsFree*P2->IsFree))
    {
    	if (!P1->IsFree) Mui = Muj;
		if (!P2->IsFree) Muj = Mui;

    	StrainRate = 2.0*vij(0)*xij(0)           , vij(0)*xij(1)+vij(1)*xij(0) , vij(0)*xij(2)+vij(2)*xij(0) ,
    				 vij(0)*xij(1)+vij(1)*xij(0) , 2.0*vij(1)*xij(1)           , vij(1)*xij(2)+vij(2)*xij(1) ,
    				 vij(0)*xij(2)+vij(2)*xij(0) , vij(1)*xij(2)+vij(2)*xij(1) , 2.0*vij(2)*xij(2)           ;
    	StrainRate = -GK * StrainRate;

    	if (VisEq==0) VI = (Mui+Muj)/(di*dj)*GK*vij;																		//Morris et al 1997
    	if (VisEq==1) VI = 4.0*(Mui+Muj)/((di+dj)*(di+dj))*GK*vij;																//Shao et al 2003
    	if (VisEq==2) VI = -(Mui+Muj)/2.0/(di*dj)* LaplaceKernel(Dimension, KernelType, rij, h)*vij;							//Real Viscosity (considering incompressible fluid)
    	if (VisEq==3) VI = -(Mui+Muj)/2.0/(di*dj)*( LaplaceKernel(Dimension, KernelType, rij, h)*vij +
    			1.0/3.0*(GK*vij + dot(vij,xij) * xij / (rij*rij) * (-GK+SecDerivativeKernel(Dimension, KernelType, rij, h) ) ) );	//Takeda et al 1994 (Real viscosity considering 1/3Mu for compressibility as per Navier Stokes but ignore volumetric viscosity)
    	if ((VisEq<0 || VisEq>3))
    	{
    	   	std::cout << "Viscosity Equation No is out of range. Please correct it and run again" << std::endl;
    		std::cout << "0 => Morris et al 1997" << std::endl;
    		std::cout << "1 => Shao et al 2003" << std::endl;
    		std::cout << "2 => Real viscosity for incompressible fluids" << std::endl;
    		std::cout << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)" << std::endl;
    	    abort();
    	}
    }
    else
    {
		if (!P1->IsFree) Mui = Muj;
		if (!P2->IsFree) Muj = Mui;
    	//Corrected velocity for fixed particle as per Morris et al 1997
    	Vec3_t vab;
    	if (P1->IsFree)
    		vab = P1->v - (2.0*P2->v-P2->vb);
   		else
   			vab = (2.0*P1->v-P1->vb) - P2->v;

    	StrainRate = 2.0*vab(0)*xij(0)           , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
				     vab(0)*xij(1)+vab(1)*xij(0) , 2.0*vab(1)*xij(1)           , vab(1)*xij(2)+vab(2)*xij(1) ,
				     vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2)           ;
		StrainRate = -GK * StrainRate;

    	if (VisEq==0) VI = (Mui+Muj)/(di*dj)*GK*vab;																		//Morris et al 1997
    	if (VisEq==1) VI = 4.0*(Mui+Muj)/((di+dj)*(di+dj))*GK*vab;																//Shao et al 2003
    	if (VisEq==2) VI = -(Mui+Muj)/2.0/(di*dj)*LaplaceKernel(Dimension, KernelType, rij, h)*vab;								//Real Viscosity (considering incompressible fluid)
    	if (VisEq==3) VI = -(Mui+Muj)/2.0/(di*dj)*( LaplaceKernel(Dimension, KernelType, rij, h)*vab +
    			1.0/3.0*(GK*vij + dot(vij,xij) * xij / (rij*rij) * (-GK+SecDerivativeKernel(Dimension, KernelType, rij, h) ) ) );	//Takeda et al 1994 (Real viscosity considering 1/3Mu for compressibility as per Navier Stokes but ignore volumetric viscosity)
    	if ((VisEq<0 || VisEq>3))
    	{
    	   	std::cout << "Viscosity Equation No is out of range. Please correct it and run again" << std::endl;
    		std::cout << "0 => Morris et al 1997" << std::endl;
    		std::cout << "1 => Shao et al 2003" << std::endl;
    		std::cout << "2 => Real viscosity for incompressible fluids" << std::endl;
    		std::cout << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)" << std::endl;
    	    abort();
    	}
    }

    // XSPH Monaghan
    if (XSPH != 0.0)
    {
        omp_set_lock(&P1->my_lock);
        P1->VXSPH		+= XSPH*mj/(0.5*(di+dj))*K*-vij;
        omp_unset_lock(&P1->my_lock);

        omp_set_lock(&P2->my_lock);
		P2->VXSPH		+= XSPH*mi/(0.5*(di+dj))*K*vij;
		omp_unset_lock(&P2->my_lock);
    }


    omp_set_lock(&P1->my_lock);
    P1->a   += -mj * ( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIij ) * GK*xij + mj*VI;
//    P1->Vis +=  mj * VI;
    if (P1->IsFree) P1->StrainRate = P1->StrainRate + mj/dj*StrainRate; else P1->StrainRate = 0.0;
    if (P1->ct==30 && Shepard)
    {
    	P1->SumDen += mj*    K;
    	P1->ZWab   += mj/dj* K;
    }
    P1->dDensity += di * (mj/dj) * dot( vij , GK*xij );
    omp_unset_lock(&P1->my_lock);


    omp_set_lock(&P2->my_lock);
    P2->a   -= -mi * ( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIji ) * GK*xij + mi*VI;
//    P2->Vis -=  mi * VI;
    if (P2->IsFree) P2->StrainRate = P2->StrainRate + mi/di*StrainRate; else P2->StrainRate = 0.0;
    if (P2->ct==30 && Shepard)
    {
    	P2->SumDen += mi*    K;
    	P2->ZWab   += mi/di* K;
    }
    P2->dDensity += dj * (mi/di) * dot( -vij , -GK*xij );
    omp_unset_lock(&P2->my_lock);
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
					LocalPairs.Push(std::make_pair(temp1, temp2));
					temp2 = Particles[temp2]->LL;
				}

				// (q1 + 1, q2 , q3)
				if (q1+1< CellNo[0])
				{
					temp2 = HOC[q1+1][q2][q3];
					while (temp2 != -1)
					{
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

	omp_set_lock(&dom_lock);
	Pairs.Push(LocalPairs);
	omp_unset_lock(&dom_lock);

}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
    {
		// Reset the pressure and the induced velocity for solid boundaries
		if (!Particles[i]->IsFree)
    	{
    		Particles[i]->vb = 0.0;
    		Particles[i]->Pressure = 0.0;
    	}

    	// The Shear-Rate, the pressure and the Bingham viscosity calculation for free particles
    	if (Particles[i]->IsFree)
    	{
            Particles[i]->Pressure = Pressure(PresEq, Cs, P0,Particles[i]->Density, Particles[i]->RefDensity);
    		Particles[i]->ShearRate = sqrt((Particles[i]->StrainRate(0,0)*Particles[i]->StrainRate(0,0) + 2.0*Particles[i]->StrainRate(0,1)*Particles[i]->StrainRate(1,0) + 2.0*Particles[i]->StrainRate(0,2)*Particles[i]->StrainRate(2,0) +
    				Particles[i]->StrainRate(1,1)*Particles[i]->StrainRate(1,1) + 2.0*Particles[i]->StrainRate(1,2)*Particles[i]->StrainRate(2,1) + Particles[i]->StrainRate(2,2)*Particles[i]->StrainRate(2,2)));
    		if (Particles[i]->T0>0.0)
			{
    			if (Particles[i]->ShearRate !=0.0)
					Particles[i]->Mu = Particles[i]->MuRef + Particles[i]->T0*(1-exp(-Particles[i]->m*Particles[i]->ShearRate))/Particles[i]->ShearRate;
				else
					Particles[i]->Mu = Particles[i]->MuRef + Particles[i]->T0*Particles[i]->m;

			}
    	}

    	//Reset to zero for all particles
    	Particles[i]->a			= a;
        Particles[i]->dDensity	= 0.0;
        Particles[i]->VXSPH		= 0.0;
        Particles[i]->ZWab		= 0.0;
        Particles[i]->StrainRate= 0.0;
        Particles[i]->SumDen	= 0.0;
        Particles[i]->Vis		= 0.0;
        Particles[i]->SumKernel	= 0.0;
    }
}

inline void Domain::PrimaryComputeAcceleration ()
{
	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (size_t k=0; k<Pairs.Size();k++)
	{
		for (size_t i=0; i<Pairs[k].Size();i++)
		{
			if (Particles[Pairs[k][i].first]->IsFree || Particles[Pairs[k][i].second]->IsFree)
			{
				if (!Particles[Pairs[k][i].first]->IsFree)
				{
					double K = Kernel(Dimension, KernelType, norm(Particles[Pairs[k][i].first]->x-Particles[Pairs[k][i].second]->x), Particles[Pairs[k][i].first]->h);
				    omp_set_lock(&Particles[Pairs[k][i].first]->my_lock);
				    Particles[Pairs[k][i].first]->SumKernel	+= K;
				    Particles[Pairs[k][i].first]->Pressure	+= Particles[Pairs[k][i].second]->Pressure * K;
				    Particles[Pairs[k][i].first]->vb		+= Particles[Pairs[k][i].second]->v * K;
				    omp_unset_lock(&Particles[Pairs[k][i].first]->my_lock);
				}
				if (!Particles[Pairs[k][i].second]->IsFree)
				{
					double K = Kernel(Dimension, KernelType, norm(Particles[Pairs[k][i].first]->x-Particles[Pairs[k][i].second]->x), Particles[Pairs[k][i].first]->h);
				    omp_set_lock(&Particles[Pairs[k][i].second]->my_lock);
				    Particles[Pairs[k][i].second]->SumKernel+= K;
				    Particles[Pairs[k][i].second]->Pressure	+= Particles[Pairs[k][i].first]->Pressure * K;
				    Particles[Pairs[k][i].second]->vb		+= Particles[Pairs[k][i].first]->v * K;
				    omp_unset_lock(&Particles[Pairs[k][i].second]->my_lock);

				}
			}
		}
	}

	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if (!Particles[i]->IsFree)
			if (Particles[i]->SumKernel!= 0.0)
			{
				Particles[i]->Pressure	= Particles[i]->Pressure/Particles[i]->SumKernel;
				Particles[i]->vb		= Particles[i]->vb/Particles[i]->SumKernel;
			}
}

inline void Domain::LastComputeAcceleration ()
{
    #pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (size_t k=0; k<Pairs.Size();k++)
	{
		for (size_t i=0; i<Pairs[k].Size();i++)
		{
			if (Particles[Pairs[k][i].first]->IsFree || Particles[Pairs[k][i].second]->IsFree)
			{
				CalcForce(Particles[Pairs[k][i].first],Particles[Pairs[k][i].second]);
			}
		}
	}
	Pairs.Clear();

	//Min time step calculation
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if (deltat > (0.25*sqrt(Particles[i]->h/norm(Particles[i]->a))) )
		{
			std::cout << "Please decrease the time step to"<< (0.25*sqrt(Particles[i]->h/norm(Particles[i]->a))) << std::endl;
			abort();
		}
}

inline void Domain::Move (double dt)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt,DomSize,TRPR,BLPF,Shepard);

}

inline void Domain::InFlowBCLeave()
{
	Array <int> DelPart;
	Array<std::pair<Vec3_t,size_t> > AddPart;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			omp_set_lock(&dom_lock);
			DelPart.Push(i);
			omp_unset_lock(&dom_lock);
		}
	if (BC.InOutFlow==2) Particles.DelItems(DelPart);

	if (BC.InOutFlow==1 || BC.InOutFlow==3)
	{
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			if(Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc1)
			{
				double temp;
				Vec3_t temp1;
				int a2,a3,b2,b3;
				temp1 = Particles[BC.InPart[i]]->x;

				if (Dimension==2) {a3=0;b3=0;}
					else {	(Particles[BC.InPart[i]]->CC[2]-1)<0 			? a3=0 				: a3=(Particles[BC.InPart[i]]->CC[2]-1);
							(Particles[BC.InPart[i]]->CC[2]+1)>=CellNo[2] 	? b3=CellNo[2]-1 	: b3=(Particles[BC.InPart[i]]->CC[2]+1);}

				(Particles[BC.InPart[i]]->CC[1]-1)<0 			? a2=0 				: a2=(Particles[BC.InPart[i]]->CC[1]-1);
				(Particles[BC.InPart[i]]->CC[1]+1)>=CellNo[1] 	? b2=CellNo[1]-1 	: b2=(Particles[BC.InPart[i]]->CC[1]+1);

				for (int q3=a3; q3<=b3; q3++)
				for (int q2=a2; q2<=b2; q2++)
				for (int q1=0; q1<Particles[BC.InPart[i]]->CC[0]; q1++)
				{
					if (HOC[q1][q2][q3]!=-1)
					{
						temp = HOC[q1][q2][q3];
						while (temp != -1)
						{
							if (Particles[temp]->IsFree && temp1(0)>Particles[temp]->x(0) && fabs(Particles[BC.InPart[i]]->x(1)-Particles[temp]->x(1))<(hmax/4.0) &&
									fabs(Particles[BC.InPart[i]]->x(2)-Particles[temp]->x(2))<(hmax/4.0) && Particles[temp]->ID==-10000) temp1=Particles[temp]->x;
							temp = Particles[temp]->LL;
						}
					}
				}
				temp1(0) -= 2.0*R;

				AddPart.Push(std::make_pair(temp1,BC.InPart[i]));
			}

		if (AddPart.Size() > DelPart.Size())
		{
			for (size_t i=0 ; i<DelPart.Size() ; i++)
			{
				Particles[DelPart[i]]->x = AddPart[i].first;
				Particles[DelPart[i]]->ID = Particles[AddPart[i].second]->ID;
				Particles[DelPart[i]]->Mass = Particles[AddPart[i].second]->Mass;
				Particles[DelPart[i]]->RefDensity = Particles[AddPart[i].second]->RefDensity;
				Particles[DelPart[i]]->h = Particles[AddPart[i].second]->h;
				Particles[AddPart[i].second]->ID = -9000;
			}
			for (size_t i=DelPart.Size() ; i<AddPart.Size() ; i++)
			{
				Particles.Push(new Particle(Particles[AddPart[i].second]->ID,AddPart[i].first,Particles[AddPart[i].second]->v,Particles[AddPart[i].second]->Mass,Particles[AddPart[i].second]->Density,Particles[AddPart[i].second]->h,false));
				Particles[AddPart[i].second]->ID = -9000;
			}
		}

		if (AddPart.Size() == DelPart.Size())
		{
			for (size_t i=0 ; i<AddPart.Size() ; i++)
			{
				Particles[DelPart[i]]->x = AddPart[i].first;
				Particles[DelPart[i]]->ID = Particles[AddPart[i].second]->ID;
				Particles[DelPart[i]]->Mass = Particles[AddPart[i].second]->Mass;
				Particles[DelPart[i]]->RefDensity = Particles[AddPart[i].second]->RefDensity;
				Particles[DelPart[i]]->h = Particles[AddPart[i].second]->h;
				Particles[AddPart[i].second]->ID = -9000;
			}
		}

		if (AddPart.Size() < DelPart.Size())
		{
			Array <int> temp1;
			for (size_t i=0 ; i<AddPart.Size() ; i++)
			{
				Particles[DelPart[i]]->x = AddPart[i].first;
				Particles[DelPart[i]]->ID = Particles[AddPart[i].second]->ID;
				Particles[DelPart[i]]->Mass = Particles[AddPart[i].second]->Mass;
				Particles[DelPart[i]]->RefDensity = Particles[AddPart[i].second]->RefDensity;
				Particles[DelPart[i]]->h = Particles[AddPart[i].second]->h;
				Particles[AddPart[i].second]->ID = -9000;
			}
			for (size_t i=AddPart.Size() ; i<DelPart.Size() ; i++)
			{
				temp1.Push(DelPart[i]);
			}
			Particles.DelItems(temp1);
		}
	}

//	if (DelPart.Size()>0 || AddPart.Size()>0) RigidParticles.Clear();
	DelPart.Clear();
	AddPart.Clear();
}

inline void Domain::InFlowBCFresh()
{
	int temp, temp1;
	int q1,q2,q3;
	if (BC.InOutFlow==1 || BC.InOutFlow==3)
	{
		BC.InPart.Clear();
		if (BC.inoutcounter==0)  BC.InFlowLoc1  = BLPF(0) + BC.cellfac*hmax;
		temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

		for (q2=0;BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0;BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=0; q1<(temp1 + 1); q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (BC.inoutcounter == 0)
						{if (Particles[temp]->IsFree && (Particles[temp]->x(0) < BC.InFlowLoc1) ) {BC.InPart.Push(temp);Particles[temp]->ID=-10000;}}
					else
						{if (Particles[temp]->IsFree && (Particles[temp]->x(0) < BC.InFlowLoc1) && (Particles[temp]->ID == -10000) ) BC.InPart.Push(temp);}

					temp = Particles[temp]->LL;
				}
			}
		}
	}

	if (BC.InOutFlow==2 || BC.InOutFlow==3)
	{
		BC.OutPart.Clear();
		if (BC.OutFlowLoc==0.0) BC.OutFlowLoc = TRPR(0) - BC.cellfac*hmax;
		temp1 = (int) (floor((BC.OutFlowLoc - BLPF(0)) / CellSize(0)));

		for (q2=0;BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0;BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=temp1;BC.Periodic[0]? (q1<(CellNo[0]-2)) : (q1<CellNo[0]) ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) > BC.OutFlowLoc) ) BC.OutPart.Push(temp);
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	// To check it is the first time that this function is running
	BC.inoutcounter = 1;

	Vec3_t vel;
	double den;

	if (BC.InPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
		{
			InCon(Particles[BC.InPart[i]]->x,vel,den,BC);
			if (norm(BC.inv)>0.0)
			{
				Particles[BC.InPart[i]]->v  = vel;
				Particles[BC.InPart[i]]->vb = vel;
			}
			if (BC.inDensity>0.0)
			{
				Particles[BC.InPart[i]]->Density  = den;
				Particles[BC.InPart[i]]->Densityb = den;
				Particles[BC.InPart[i]]->RefDensity = BC.inDensity;
			}
		}

	double temp11;
	if (BC.MassConservation) temp11= BC.InPart.Size()*BC.inv(0)/BC.OutPart.Size();

	if (BC.OutPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.OutPart.Size() ; i++)
		{
			OutCon(Particles[BC.OutPart[i]]->x,vel,den,BC);
			if (temp11<vel(0) && BC.MassConservation) vel(0) = temp11;
			if (norm(BC.outv)>0.0)
			{
				Particles[BC.OutPart[i]]->v  = vel;
				Particles[BC.OutPart[i]]->vb = vel;
			}
			if (BC.outDensity>0.0)
			{
				Particles[BC.OutPart[i]]->Density  = den;
				Particles[BC.OutPart[i]]->Densityb = den;
				Particles[BC.OutPart[i]]->RefDensity = BC.outDensity;
			}
		}

}

inline void Domain::InFlowBCReset()
{
	if (BC.InPart.Size()>0)
		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
		{
				Particles[BC.InPart[i]]->a = 0.0;
				Particles[BC.InPart[i]]->dDensity  = 0.0;
				Particles[BC.InPart[i]]->SumDen  = 0.0;
				Particles[BC.InPart[i]]->ZWab  = 0.0;
		}

	if (BC.OutPart.Size()>0)
		#pragma omp parallel for schedule(static) num_threads(Nproc)
		for (size_t i=0 ; i<BC.OutPart.Size() ; i++)
		{
				Particles[BC.OutPart[i]]->a = 0.0;
				Particles[BC.OutPart[i]]->dDensity  = 0.0;
				Particles[BC.OutPart[i]]->SumDen  = 0.0;
				Particles[BC.OutPart[i]]->ZWab  = 0.0;
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
				Particles[i]->vb = vel;
    		}
    		if (Particles[i]->IsFree && BC.allDensity>0.0)
    		{
				Particles[i]->Density = den;
				Particles[i]->Densityb = den;
				Particles[i]->RefDensity = BC.allDensity;
    		}
    	}
    }
}

inline void Domain::InitialChecks()
{
    if (KernelType==4) Cellfac = 3.0; else Cellfac = 2.0;

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

    if (BCDensityUpdate==false)
    {
		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0 ; i<Particles.Size() ; i++)
			if (!Particles[i]->IsFree)
				Particles[i]->DensityUpdate = false;
    }


    InitialChecks();
    CellInitiate();
    ListGenerate();
    PrintInput(TheFileKey);
    WholeVelocity();

    while (Time<tf && idx_out<=maxidx)
    {
    	StartAcceleration(Gravity);

    	InFlowBCFresh();

    	MainNeighbourSearch();

    	PrimaryComputeAcceleration();

    	GeneralBefore(*this);

    	LastComputeAcceleration();

    	GeneralAfter(*this);

    	InFlowBCReset();

    	Move(dt);

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

       if (BC.InOutFlow>0 && !BC.Periodic[0]) InFlowBCLeave(); else CheckParticleLeave ();
//       if (!(BC.Periodic[0] || BC.Periodic[1] || BC.Periodic[2])) CheckParticleLeave ();

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

    oss << "Viscosity Equation = ";
    if (Alpha!=0.0 || Beta!=0.0) oss << "Artificial Viscosity by Monaghan\nAlpha = "<<Alpha<<"   Beta = "<<Beta <<"\n";
    else
	{
    	switch (VisEq)
		{
			case 0:
				oss << "0 => Morris et al 1997\n";
				break;
			case 1:
				oss << "1 => Shao et al 2003\n";
				break;
			case 2:
				oss << "2 => Real viscosity for incompressible fluids\n";
				break;
			case 3:
				oss << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)\n";
				break;
		}
	}

    oss << "Equation of State = ";
	switch (PresEq)
    {
		case 0:
			oss << "0 => P0+(Cs*Cs)*(Density-Density0)\n";
			break;
		case 1:
			oss << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)\n";
			break;
		case 2:
			oss << "2 => (Cs*Cs)*Density\n";
			break;
    }
	oss << "Cs = "<<Cs<<" m/s\n";
	oss << "P0 = "<<P0<<" Pa\n";
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
    t1 = 0.25*hmax/(Cs);
    t2 = 0.125*hmax*hmax*rhomax/MuMax;

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
    oss << "No-Slip Boundary Condition = " << (NoSlip ? "True" : "False")  << "\n";
    oss << "Shepard Filter for Density = " << (Shepard ? "True" : "False") << "\n";
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
    dims[0] = 6;



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
