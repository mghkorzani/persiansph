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

#include "Domain.h"

using std::cout;
using std::endl;
void NewUserOutput(SPH::Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = Particles->ZWab;
}

int main(int argc, char **argv) try
{
        SPH::Domain	dom;

        dom.Dimension	= 2;
        dom.Nproc	= 24;
    	dom.KernelType	= 4;
    	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs;
    	double H,L,n;

    	H	= 0.02;
    	L	= 0.2;
    	n	= 20.0;
    	rho	= 1000.0;
    	K	= 3.25e6;
    	G	= 7.15e5;
    	dx	= H / n;
    	h	= dx*1.1; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.2*h/(Cs));
        cout<<timestep<<endl;
        cout<<Cs<<endl;
        dom.DomMax(1) = 0.10;
        dom.DomMin(1) =-0.10;
        dom.DomMax(0) = 0.22;


     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , -H/2.0-4.0*dx , 0.0 ), 1.5*L + dx/10.0 , H + 8.0*dx + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

     	double Vy, Vf, k, x, y;
     	Vf = 0.01;
     	k = 1.875/0.2;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
    		y = dom.Particles[a]->x(1);
    		if (x>0.0 && y>H/2.0) dom.Particles[a]->ID=5;
    		if (x>0.0 && y<-H/2.0) dom.Particles[a]->ID=5;
    		if (x<=0.0 && y>H/2.0)
    		{
    			dom.Particles[a]->ID=2;
    			dom.Particles[a]->IsFree=false;
    			dom.Particles[a]->NoSlip= true;
    		}
    		if (x<=0.0 && y<-H/2.0)
    		{
     			dom.Particles[a]->ID=2;
    			dom.Particles[a]->IsFree=false;
    			dom.Particles[a]->NoSlip= true;
    		}
    	}
    	dom.DelParticles(5);

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G		= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Alpha		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;

    		x = dom.Particles[a]->x(0);
    		if (x>0.0)
    		{
    			Vy = Cs*Vf*((cos(k*L)+cosh(k*L))*(cosh(k*x)-cos(k*x))+(sin(k*L)-sinh(k*L))*(sinh(k*x)-sin(k*x)))/((cos(k*L)+cosh(k*L))*(cosh(k*L)-cos(k*L))+(sin(k*L)-sinh(k*L))*(sinh(k*L)-sin(k*L)));
    			dom.Particles[a]->v(1)  = Vy;
    			dom.Particles[a]->vb(1) = Vy;
    		}
    	}
	double T, Nu;
	Nu	= (3.0*K-2.0*G)/(6.0*K+2.0*G);
	T	= std::sqrt((4.0*M_PI*M_PI*6.0*rho*(1-Nu))/(G*H*H*k*k*k*k));
	cout<<"Nu = "<<Nu<<endl;
	cout<<"T  = "<<T<<endl;
	cout<<"TCs/L  = "<<T*Cs/L<<endl;

	dom.OutputName[0]	= "ZWab";
        dom.UserOutput		= & NewUserOutput;
	


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.002,"test06",999);
        return 0;
}
MECHSYS_CATCH
