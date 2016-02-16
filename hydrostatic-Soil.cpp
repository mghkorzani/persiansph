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

double check = 0.0;
double DampF,DampS;
double t = 0.2;
void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<t)
    #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 1) domi.Particles[i]->a -= DampF * domi.Particles[i]->v;
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;
	}
}

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.KernelType	= 0;
	dom.Nproc		= 24;

//	dom.Shepard		= true;
//	dom.XSPH		= 0.5;
	dom.Scheme		= 0;
	dom.Gravity		= 0.0,-9.81,0.0;
//	dom.TI			= 0.3;
	dom.GeneralAfter= & UserDamping;
	dom.BC.Periodic[0] = true;

	double xb,yb,h,rhoF,rhoS,CsF,CsS;
	double dx,HF,HS,L,K,G,Nu,E,T1,T2,T;

	HF				= 1.4;
	HS				= 1.0;
	L				= 0.4;
	rhoF			= 998.21;
	rhoS			= 2038.7;
	dx				= 0.02;
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
   	DampS 			= 0.04*sqrt(E/(rhoS*h*h));
	DampF 			= 0.04*CsF/h;

	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , -3.0*dx , 0.0 ), L + dx/10.0 , HF + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Alpha		= 1.0;
		dom.Particles[a]->Beta		= 1.0;
		dom.Particles[a]->Material	= 1;
		dom.Particles[a]->Cs		= CsF;
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
		}
		if (dom.Particles[a]->ID==1) dom.Particles[a]->Density  = rhoF*((1+9.81*(HF-dom.Particles[a]->x(1))/(CsF*CsF)));

	}


	dom.AddBoxLength(3 ,Vec3_t ( -L/2.0 , -3.0*dx , 0.0 ), L + dx/10.0 , HS + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoS, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==3)
		{
			dom.Particles[a]->n			= 0.000001;
			dom.Particles[a]->k			= 1.0;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G			= G;
			dom.Particles[a]->K			= K;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Fail		= 0;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			dom.Particles[a]->RhoF		= rhoF;
			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);
			if (yb<0.0)
			{
				dom.Particles[a]->ID=4;
				dom.Particles[a]->IsFree=false;
				dom.Particles[a]->NoSlip=true;
			}
		}
	}

	dom.Solve(/*tf*/50000.0,/*dt*/T,/*dtOut*/0.01,"test06",10000);
	return 0;
}
MECHSYS_CATCH
