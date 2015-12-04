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

#include "./Source/Domain.h"
#include "./Source/Interaction.h"


using std::cout;
using std::endl;
using std::ifstream;


int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.KernelType		= 0;
	dom.Nproc			= 24;
	dom.Scheme			= 0;
	dom.VisEq			= 0;
	dom.Gravity			= 0.0,-9.81,0.0;

	double h,rhoF,rhoS,CsF,CsS;
	double dx,HF,K,G,Nu,E,T1,T2,T;

	HF				= 1;
	rhoF			= 998.21;
	rhoS			= 2038.7;
	dx				= 0.004;
	h				= dx*1.3;
	E				= 15.0e6;
	Nu				= 0.3;
	K				= E/(3.0*(1.0-2.0*Nu));
	G				= E/(2.0*(1.0+Nu));
    CsS				= 600.0/6.0; //Very important
	CsF				= 10.0*sqrt(2.0*9.81*HF);
	dom.InitialDist	= dx;
	T1				= (0.1*h/CsF);
    T2				= (0.2*h/CsS);
    T				= std::min(T1,T2);
    std::cout<<"Tf = "<<T1<<std::endl;
    std::cout<<"Ts = "<<T2<<std::endl;
    std::cout<<"T  = "<<T<<std::endl;

	dom.AddBoxLength(1 ,Vec3_t ( 0.0           , -3.0*dx , 0.0 ), 0.1 + 6.0*dx + dx/10.0 , 1.0 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);
	dom.AddBoxLength(2 ,Vec3_t ( 0.1 + 6.0*dx  , -3.0*dx , 0.0 ), 0.2 + 6.0*dx + dx/10.0 , 0.5 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);
	dom.AddBoxLength(4 ,Vec3_t ( 0.1 + 6.0*dx  , 0.2     , 0.0 ), 0.2 + 6.0*dx + dx/10.0 , 0.2 + dx/10.0           ,  0 , dx/2.0 ,rhoS, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==1) dom.Particles[a]->Density  = rhoF*((1+9.81*(1.0-dom.Particles[a]->x(1))/(CsF*CsF)));
		if (dom.Particles[a]->ID==2) dom.Particles[a]->Density  = rhoF*((1+9.81*(0.5-dom.Particles[a]->x(1))/(CsF*CsF)));
		if (dom.Particles[a]->ID == 1 || dom.Particles[a]->ID == 2)
		{
			dom.Particles[a]->Alpha		= 0.05;
			dom.Particles[a]->Mu		= 1.002e-3;
			dom.Particles[a]->MuRef		= 1.002e-3;
			dom.Particles[a]->Material	= 1;
			dom.Particles[a]->Cs		= CsF;
			if (dom.Particles[a]->x(0)<3.0*dx)
			{
				dom.Particles[a]->ID		= 3;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
			}
			if (dom.Particles[a]->x(1)<0.0)
			{
				dom.Particles[a]->ID		= 3;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
			}
			if (dom.Particles[a]->x(0)>(0.1+3.0*dx) && dom.Particles[a]->x(0)<(0.1+9.0*dx))
			{
				dom.Particles[a]->ID		= 3;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
			}
			if (dom.Particles[a]->x(0)>(0.1+3.0*dx) && dom.Particles[a]->x(0)<(0.1+9.0*dx))
			{
				dom.Particles[a]->ID		= 3;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
			}
			if (dom.Particles[a]->x(0)>(0.3+9.0*dx) && dom.Particles[a]->x(0)<(0.3+12.0*dx))
			{
				dom.Particles[a]->ID		= 3;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
			}
		}
		if (dom.Particles[a]->x(0)>(0.1+3.0*dx) && dom.Particles[a]->x(0)<(0.1+9.0*dx)
				&& dom.Particles[a]->x(1)>0.0  && dom.Particles[a]->x(1)<0.15)
		{
			dom.Particles[a]->ID		= 5;
			dom.Particles[a]->IsFree	= true;
		}
		if (dom.Particles[a]->ID == 4)
		{
			dom.Particles[a]->n			= 0.09;
			dom.Particles[a]->k			= 5.0e-3;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G			= G;
			dom.Particles[a]->K			= K;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Fail		= 0;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			dom.Particles[a]->IsFree	= false;
		}
	}



	dom.Solve(/*tf*/50000.0,/*dt*/T,/*dtOut*/0.01,"test06",10000);
	return 0;
}
MECHSYS_CATCH
