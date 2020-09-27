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

void NewUserOutput(SPH::Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = Particles->ZWab;
}


using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char **argv) try
{
        SPH::Domain	dom;

        dom.Dimension	= 2;
        dom.Nproc	= 24;
    	dom.VisEq	= 0;
    	dom.KernelType	= 0;
    	dom.Scheme	= 0;
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;
	dom.XSPH	= 0.5;


    	double	h,dx,g;
    	dx	= 0.001;
    	h	= dx*1.4;
	g	= norm(dom.Gravity);

	double	L,L1,H,H1,RhoF,CsW,MuW,t1,Tr;
	L	= 0.25;
	L1	= 0.1;
	H	= 0.14;
	H1	= 0.079+2.0*dx;
	Tr	= 0.005;
	RhoF	= 1000.0;
	CsW	= 10.0*sqrt(2.0*g*H)*2.0;
	MuW	= 1.3e-3;
        t1	= (0.2*h/(CsW));
	cout<<"CsW = "<<CsW<<endl;
	cout<<"t1  = "<<t1<<endl;



    	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , 0.0 ), L + 3.0*dx + dx/10.0 , H + 6.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

    	double yb,xb;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		xb				= dom.Particles[a]->x(0);
    		yb				= dom.Particles[a]->x(1);

    		dom.Particles[a]->Cs		= CsW;
    		dom.Particles[a]->Alpha		= 0.1;
    		dom.Particles[a]->Beta		= 0.1;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Mu		= MuW;
    		dom.Particles[a]->MuRef		= MuW;
    		dom.Particles[a]->Material	= 1;
//    		dom.Particles[a]->Shepard	= true;

    		if (yb<0.0 || xb<0.0)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    		if (xb>L1 && xb<=(L1+3.0*dx) && yb>H1)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    		if (xb>L1 && dom.Particles[a]->ID == 1)
    			dom.Particles[a]->ID		= 3;

    		if (yb>H && dom.Particles[a]->ID == 1)
    			dom.Particles[a]->ID		= 3;

//    		if (dom.Particles[a]->ID == 1)
//		{
//	    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H-yb)/(CsW*CsW)),(1.0/7.0));
//  			dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H-yb)/(CsW*CsW)),(1.0/7.0));
//		}
   	}
	
	dom.DelParticles(3);



	//Rubber gate definition
        double	RhoR,K,G,CsR,t2;
    	RhoR	= 1100.0;
    	K	= 2.00e7;
    	G	= 4.27e6;
        CsR	= sqrt(K/RhoR);
//	CsR	= 100.0;
        t2	= (0.2*h/(CsR));
	cout<<"CsR = "<<CsR<<endl;
	cout<<"t2  = "<<t2<<endl;

     	dom.AddBoxLength(3 ,Vec3_t ( L1 , 0.0 , 0.0 ), Tr + 3.0*dx + dx/10.0 , H1 + 3.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoR, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
		if (dom.Particles[a]->ID == 3)
		{
	    		xb = dom.Particles[a]->x(0);
	    		yb = dom.Particles[a]->x(1);

	    		dom.Particles[a]->G		= G;
	    		dom.Particles[a]->PresEq	= 0;
	    		dom.Particles[a]->Cs		= CsR;
	    		dom.Particles[a]->Material	= 2;
	    		dom.Particles[a]->Shepard	= false;
	    		dom.Particles[a]->Alpha		= 1.0;
	    		dom.Particles[a]->TI		= 0.3;
	    		dom.Particles[a]->TIInitDist	= dx;

	    		if (yb>H1)
	    		{
	    			dom.Particles[a]->ID	= 4;
	    			dom.Particles[a]->IsFree= false;
	    			dom.Particles[a]->NoSlip= true;
	    		}
	    		if (xb>(L1+Tr) && yb<=H1)
	      			dom.Particles[a]->ID	= 5;

	    		if (yb<2.0*dx)
	      			dom.Particles[a]->ID	= 5;
		}
    	}
    	dom.DelParticles(5);

	double	t;
	t	= std::min(t1,t2);
	cout<<"t   = "<<t<<endl;

	dom.OutputName[0]	= "ZWab";
        dom.UserOutput		= & NewUserOutput;

//	dom.WriteXDMF("maz");
   	dom.Solve(/*tf*/0.5,/*dt*/t,/*dtOut*/0.001,"test",1000);

        return 0;
}
MECHSYS_CATCH
