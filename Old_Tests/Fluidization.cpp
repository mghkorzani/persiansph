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
#include "Interaction.h"
	double g,dx,RhoF,CsW,DampTime,DampF,DampS,L,H,Q,V,q1=0.0,Muw,d;
	int	check=0;
void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = V,0.0,0.0;
	Den = RhoF*pow((1+7.0*(g*(H-4.25*dx)+0.5*V*V)/(CsW*CsW)),(1.0/7.0));
}
void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = V,0.0,0.0;
	Den = RhoF*pow((1+7.0*(g*(4.25*dx)+0.5*V*V)/(CsW*CsW)),(1.0/7.0));
}

void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<=DampTime)
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
			if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 1) domi.Particles[i]->a  = 0.0;
			if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;
		}
	}
	else
	{
		if (Q<=0.12)
		{
			Q	= 0.01*(domi.Time-DampTime);
			V	= (Q/1000.0)/(0.15*0.15);
		}
		if (check==0)
		{
			#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
			for (size_t i=0; i<domi.Particles.Size(); i++)
			{
				if (domi.Particles[i]->IsFree)
				{ 
					domi.Particles[i]->v = 0.0;
					domi.Particles[i]->vb = 0.0;
				}
			}
			domi.SWIType		= 0;
			domi.BC.InOutFlow	= 3;
			domi.BC.inv		= V,0.0,0.0;
			domi.BC.outv		= V,0.0,0.0;
			domi.BC.inDensity	= 998.21;
			domi.BC.outDensity	= 998.21;
			domi.InCon		= & UserInFlowCon;
			domi.OutCon		= & UserOutFlowCon;
		}
		check = 1;

		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
			if (domi.Particles[i]->Material == 1)
			{ 
				domi.Particles[i]->a(1) = 0.0;
			}
		}

		if (Q>q1)
		{
			std::cout<<"Q = "<<Q<<std::endl;
			std::cout<<"Re = "<<(d*RhoF*V/Muw)<<std::endl;
			q1 += 0.01;
		}
	}
		
}

void UserDamping1(SPH::Domain & domi)
{
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3 && !(domi.Particles[i]->IsSat)) 
		{
			domi.Particles[i]->Mass		= domi.Particles[i]->V*(domi.Particles[i]->RefDensity - domi.Particles[i]->RhoF);
			domi.Particles[i]->Density	= domi.Particles[i]->Density - domi.Particles[i]->RhoF;
			domi.Particles[i]->Densityb	= domi.Particles[i]->Densityb - domi.Particles[i]->RhoF;
			domi.Particles[i]->RefDensity	= domi.Particles[i]->RefDensity - domi.Particles[i]->RhoF;
			domi.Particles[i]->IsSat	= true;
			domi.Particles[i]->SatCheck	= true;
		}
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3 ) 
		{
			domi.Particles[i]->IsSat	= true;
			domi.Particles[i]->SatCheck	= true;
		}
	}
}


void NewUserOutput(SPH::Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = Particles->ZWab;
	Prop2 = Particles->n;
	Prop3 = Particles->k;
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
	dom.SWIType	= 3;
    	dom.Scheme	= 0;
    	dom.Gravity	= -9.81, 0.0 , 0.0 ;
	g		= norm(dom.Gravity);

        dom.BC.Periodic[1]	= true;

        dom.GeneralAfter	= & UserDamping;
        dom.GeneralBefore	= & UserDamping1;


    	double h,t,t1,t2,h1,h2,h3,h0;
    	dx	= 0.001;
    	h	= dx*1.2;
	dom.InitialDist	= dx;
	dom.BC.cellfac	= 4.25;

//	L	= 0.15;
	L	= 0.06;
	h0	= 0.0;
	h1	= 0.0125;
	h2	= 0.0125;
	h3	= 0.025;

	H	= h0+h1+h2+h3+0.025;
	Q	= 0.05;
	V	= (Q/1000.0)/(0.15*0.15);
	RhoF	= 1000.0;
	CsW	= 10.0*sqrt(2.0*g*H);
	Muw	= 1.0e-3;
        t1	= (0.25*h/(CsW));



    	dom.AddBoxLength(1 ,Vec3_t ( -6.1*dx , -L/2.0 , 0.0 ), H + dx/10.0 , L + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

    	double yb,xb;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		xb=dom.Particles[a]->x(0);
    		yb=dom.Particles[a]->x(1);

    		dom.Particles[a]->Cs		= CsW;
    		dom.Particles[a]->Alpha		= 0.05;
//    		dom.Particles[a]->Beta		= 0.05;
    		dom.Particles[a]->PresEq	= 1;
    		dom.Particles[a]->Mu		= Muw;
    		dom.Particles[a]->MuRef		= Muw;
    		dom.Particles[a]->Material	= 1;
    		dom.Particles[a]->Shepard	= true;
    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H-xb)/(CsW*CsW)),(1.0/7.0));
		dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H-xb)/(CsW*CsW)),(1.0/7.0));
/*
    		if (xb<-L/2.0 || xb>L/2.0)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
//    			dom.Particles[a]->NoSlip	= true;
    		}
*/
   	}

	double Nu,E,K,G,CsS,RhoS,c,Phi,Psi,n1,n2,n3,d1,d2,d3;

	Nu	= 0.3;
	E	= 5.0e6;
	K	= E/(3.0*(1.0-2.0*Nu));
	G	= E/(2.0*(1.0+Nu));
	n1	= 0.36;
	n2	= 0.36;
	n3	= 0.46;
	RhoS	= 2500.0*(1.0-n3)+n3*RhoF;
	CsS	= sqrt(K/(RhoS-RhoF));
	c	= 0.0;
	Phi	= 22.0;
	Psi	= 0.0;
	d1	= 0.003;
	d2	= 0.002;
	d3	= 0.0005125;
	d	= d3;
        t2	= (0.25*h/(CsS));

        std::cout<<"t2 = "<<t2<<std::endl;

        std::cout<<"CsS  = "<<CsS<<std::endl;
        std::cout<<"RhoS = "<<RhoS<<std::endl;
        std::cout<<"Phi  = "<<Phi<<std::endl;
        std::cout<<"C    = "<<c<<std::endl;
        std::cout<<"Q    = "<<Q<<std::endl;
        std::cout<<"V    = "<<V<<std::endl;

	dom.AddBoxLength(5 ,Vec3_t ( h0 , -L/2.0 , 0.0 ), h0+h1+h2+h3 + dx/10.0 , L + dx/10.0 ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==5)
		{
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			if (c>0.0)
			{
				dom.Particles[a]->TI	= 0.5;
				dom.Particles[a]->TIn	= 2.55;
			}
			dom.Particles[a]->TIInitDist	= dx;
			dom.Particles[a]->Shepard	= true;
			dom.Particles[a]->VarPorosity	= true;
			dom.Particles[a]->SeepageType	= 1;	
			dom.Particles[a]->RhoF		= RhoF;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G		= G;
			dom.Particles[a]->K		= K;
			dom.Particles[a]->Fail		= 3;	//non-associated flow rule 
			dom.Particles[a]->c		= c;
			dom.Particles[a]->phi		= Phi/180.0*M_PI;
			dom.Particles[a]->psi		= Psi/180.0*M_PI;
	
			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);
/*
	    		if (xb<-L/2.0 || xb>L/2.0)
	    		{
	    			dom.Particles[a]->ID		= 8;
				dom.Particles[a]->d		= d1*1.0e20;
				dom.Particles[a]->n0		= 1.0;
				dom.Particles[a]->n		= 1.0;
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
	    		}
*/

	    		if (xb>(h0) && xb<=(h0+h1) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 5;
				dom.Particles[a]->d		= d1;
				dom.Particles[a]->n0		= n1;
				dom.Particles[a]->n		= n1;
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
	    		}
	    		if (xb>(h0+h1) && xb<=(h0+h1+h2) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 6;
				dom.Particles[a]->d		= d2;
				dom.Particles[a]->n0		= n2;
				dom.Particles[a]->n		= n2;
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
	    		}
	    		if (xb>(h0+h1+h2) && xb<=(h0+h1+h2+h3) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 7;
				dom.Particles[a]->d		= d3;
				dom.Particles[a]->n0		= n3;
				dom.Particles[a]->n		= n3;
	    		}
		}
	}
//	dom.DelParticles(10);

   	DampF	= 0.02*CsW/h;
  	DampS	= 0.02*sqrt(E/(RhoS*h*h));
    	DampTime= 0.02;
//    	DampTime= 0.0;

        t	= std::min(t1,t2);
        t	= 2.0e-6;
        std::cout<<"t1 = "<<t1<<std::endl;
        std::cout<<"t2 = "<<t2<<std::endl;
        std::cout<<"t  = "<<t<<std::endl;

	dom.OutputName[0]	= "ZWab";
	dom.OutputName[1]	= "Porosity";
	dom.OutputName[2]	= "Permeability";
        dom.UserOutput		= & NewUserOutput;

   	dom.Solve(/*tf*/99.8,/*dt*/t,/*dtOut*/0.1,"test",100000);
        return 0;
}
MECHSYS_CATCH
