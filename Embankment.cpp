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
double t = 0.5;
void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<t)
    #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;
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
//	dom.GeneralAfter= & UserDamping;
//	dom.BC.Periodic[0] = true;

	double xb,yb,h,dx,T;

	dx		= 0.15;
	h		= dx*1.3;
	dom.InitialDist	= dx;

	double rhoF,Mu,CsF,Tf;
	rhoF	= 998.21;
	Mu		= 1.002e-3;
	CsF		= 10.0*sqrt(2.0*9.81*4.5);
	Tf		= (0.2*h/CsF);
	DampF 	= 0.04*CsF/h;
	dom.VisEq= 0;


	dom.AddBoxLength(1 ,Vec3_t ( -16.0 - 3.0*dx , -3.0*dx , 0.0 ), 32.0 + 6.0*dx + dx/10.0 , 4.5 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Alpha		= 0.03;
		dom.Particles[a]->Mu		= Mu;
		dom.Particles[a]->MuRef		= Mu;
		dom.Particles[a]->Material	= 1;
		dom.Particles[a]->Cs		= CsF;
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID	= 2;
			dom.Particles[a]->NoSlip= true;
			dom.Particles[a]->IsFree= false;
		}
		if (xb<-16.0)
		{
			dom.Particles[a]->ID	= 2;
			dom.Particles[a]->NoSlip= true;
			dom.Particles[a]->IsFree= false;
		}
		if (yb>(-0.5*(xb-11.0)) && dom.Particles[a]->ID == 1)
			dom.Particles[a]->ID	= 5;

		if (dom.Particles[a]->ID==1 && xb<=2.0) dom.Particles[a]->Density  = rhoF*((1+9.81*(4.5-dom.Particles[a]->x(1))/(CsF*CsF)));
		if (dom.Particles[a]->ID==1 && xb>2.0 ) dom.Particles[a]->Density  = rhoF*((1+9.81*((-0.5*(xb-11.0))-dom.Particles[a]->x(1))/(CsF*CsF)));
	}

	double K,G,Nu,E,rhoS,CsS,Phi,c,Psi,k,Ts;

	rhoS	= 2038.7;
	Phi		= 25.0;
	Psi		= 0.0;
	k		= 2.0;
	c		= 10.0e3;
	E		= 25.0e6;
	Nu		= 0.3;
	K		= E/(3.0*(1.0-2.0*Nu));
	G		= E/(2.0*(1.0+Nu));
    CsS		= 600.0/6.0; //Very important
    Ts		= (0.2*h/CsS);
   	DampS	= 0.05*sqrt(E/(rhoS*h*h));

	dom.AddBoxLength(3 ,Vec3_t ( -16.0 , -3.0*dx , 0.0 ), 32.0 + dx/10.0 , 5.0 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoS, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==3)
		{
			dom.Particles[a]->n			= 1.0;
			dom.Particles[a]->k			= k;
			dom.Particles[a]->RhoF		= rhoF;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G			= G;
			dom.Particles[a]->K			= K;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Fail		= 3;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			dom.Particles[a]->TI		= 0.5;
			dom.Particles[a]->TIn		= 2.55;
    		dom.Particles[a]->c			= c;
    		dom.Particles[a]->phi		= Phi/180.0*M_PI;
    		dom.Particles[a]->psi		= Psi/180.0*M_PI;
			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);
			if (yb<0.0)
			{
				dom.Particles[a]->ID	= 4;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= true;
			}
			if (yb>(-0.5*(xb-11.0)) && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 5;
			if (yb>( 0.5*(xb+11.0)) && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 5;
		}
	}
	dom.DelParticles(5);

    T		= std::min(Tf,Ts);
    std::cout<<"Tf = "<<Tf<<std::endl;
    std::cout<<"Ts = "<<Ts<<std::endl;
    std::cout<<"T  = "<<T<<std::endl;

	dom.Solve(/*tf*/50000.0,/*dt*/T,/*dtOut*/0.1,"test06",10000);
	return 0;
}
MECHSYS_CATCH
