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
	double g,dx,RhoF,CsW,DampTime,DampF,DampS,L,H,Hs;
	int	check=0;

void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<DampTime)
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
			if (domi.Particles[i]->IsFree && domi.Particles[i]->ID == 4)		domi.Particles[i]->a = 0.0;
			if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 1)	domi.Particles[i]->a -= DampF * domi.Particles[i]->v;
			if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3)	domi.Particles[i]->a -= DampS * domi.Particles[i]->v;
		}
	}
	else
	{
		if (check==0)
		{
			domi.SWIType	= 1;
			#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
			for (size_t i=0; i<domi.Particles.Size(); i++)
			{
				if (domi.Particles[i]->IsFree)
				{ 
//					if (domi.Particles[i]->Material == 1)	domi.Particles[i]->LES	= true;
					if (domi.Particles[i]->ID == 4)	domi.Particles[i]->ID == 1;
					domi.Particles[i]->v = 0.0;;
					domi.Particles[i]->vb = 0.0;;
				}
			}
			domi.DelParticles(7);
			domi.DelParticles(3);

		}
		check = 1;
	}
		
}


void NewUserOutput(SPH::Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = Particles->ZWab;
	Prop2 = Particles->n;
	Prop3 = Particles->Mu - Particles->MuRef;
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
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;
	g		= norm(dom.Gravity);

        dom.GeneralAfter	= & UserDamping;


    	double h,t,t1,t2,Muw;
    	dx	= 0.005;
    	h	= dx*1.2;
	dom.InitialDist	= dx;

	L	= 6.0;
	H	= 0.35;
	Hs	= 0.10;

	RhoF	= 1000.0;
	CsW	= 10.0*sqrt(2.0*g*(H+Hs))*2.0;
	Muw	= 1.0e-3;
        t1	= (0.25*h/(CsW));



    	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-4.0*dx , -4.0*dx , 0.0 ), L + 8.0*dx + dx/10.0 , H + Hs + 8.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

    	double yb,xb;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		xb=dom.Particles[a]->x(0);
    		yb=dom.Particles[a]->x(1);

    		dom.Particles[a]->Cs		= CsW;
    		dom.Particles[a]->Alpha		= 0.05;
    		dom.Particles[a]->Beta		= 0.05;
    		dom.Particles[a]->PresEq	= 1;
    		dom.Particles[a]->Mu		= Muw;
    		dom.Particles[a]->MuRef		= Muw;
    		dom.Particles[a]->Material	= 1;
//    		dom.Particles[a]->LES		= true;
//    		dom.Particles[a]->Shepard	= true;



     		if (xb<(-L/2.0) || xb>(L/2.0))
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}

	    	if (xb>0.0 && xb<4.0*dx && yb>Hs && dom.Particles[a]->ID==1)
		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}

   		if (xb>0.0 && yb>Hs && dom.Particles[a]->ID==1)
    			dom.Particles[a]->ID		= 10;
   		if (yb>(H+Hs) && dom.Particles[a]->ID==1)
    			dom.Particles[a]->ID		= 10;


    		if (dom.Particles[a]->ID==1 && xb<0.0)
  		{
	    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H+Hs-yb)/(CsW*CsW)),(1.0/7.0));
    			dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H+Hs-yb)/(CsW*CsW)),(1.0/7.0));
    		}
    		if (dom.Particles[a]->ID==1 && xb>0.0)
  		{
	    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(Hs-yb)/(CsW*CsW)),(1.0/7.0));
    			dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(Hs-yb)/(CsW*CsW)),(1.0/7.0));
    		}

	    	if (xb>0.0 && xb<4.0*dx && yb<Hs && dom.Particles[a]->ID==1)
		{
    			dom.Particles[a]->ID		= 4;
    		}


   	}
	dom.DelParticles(10);

	double Nu,E,K,G,CsS,RhoS,c,Phi,Psi,n,d;

	Nu	= 0.3;
	E	= 5.0e6;
	K	= E/(3.0*(1.0-2.0*Nu));
	G	= E/(2.0*(1.0+Nu));
	n	= 0.42;
	RhoS	= 1580.0*(1.0-n)+n*RhoF;
	CsS	= sqrt(K/(RhoS-RhoF));
	c	= 0.0;
	Phi	= 38.0;
	Psi	= 0.0;
	d	= 0.00392;
        t2	= (0.05*h/CsS);

        std::cout<<"CsS  = "<<CsS<<std::endl;
        std::cout<<"RhoS = "<<RhoS<<std::endl;
        std::cout<<"Phi  = "<<Phi<<std::endl;
        std::cout<<"C    = "<<c<<std::endl;

	dom.AddBoxLength(5 ,Vec3_t ( -L/2.0-4.0*dx , -4.0*dx , 0.0 ), L + 8.0*dx + dx/10.0 , H + Hs + 8.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, false );

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
			dom.Particles[a]->d		= d;
//	    		dom.Particles[a]->Shepard	= true;
//			dom.Particles[a]->VarPorosity	= true;
			dom.Particles[a]->SeepageType	= 1;	
			dom.Particles[a]->n0		= n;
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

			if (xb<(-L/2.0) || xb>(L/2.0))
			{
				dom.Particles[a]->ID		= 6;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
				dom.Particles[a]->d		= d*1.0e20;
			}
			if (yb<0.0)
			{
				dom.Particles[a]->ID		= 6;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
				dom.Particles[a]->d		= d*1.0e20;
			}

	    		if (xb>0.0 && xb<4.0*dx && yb>Hs && yb<(Hs+0.1+4.0*dx) && dom.Particles[a]->ID==5)
			{
				dom.Particles[a]->ID		= 7;
				dom.Particles[a]->IsFree	= false;
				dom.Particles[a]->NoSlip	= true;
				dom.Particles[a]->d		= d*1.0e20;
			}

	    		if (((xb<0.0 && yb>(Hs+0.1)) || (xb>0.0 && yb>Hs))  && dom.Particles[a]->ID==5)
	    			dom.Particles[a]->ID		= 10;
		}
	}
	dom.DelParticles(10);

   	DampF	= 0.02*CsW/h;
  	DampS	= 0.02*sqrt(E/(RhoS*h*h));
    	DampTime= 0.1;

        t	= std::min(t1,t2);
        std::cout<<"t1 = "<<t1<<std::endl;
        std::cout<<"t2 = "<<t2<<std::endl;
        std::cout<<"t  = "<<t<<std::endl;

	dom.OutputName[0]	= "ZWab";
	dom.OutputName[1]	= "Porosity";
	dom.OutputName[2]	= "MuLES";
        dom.UserOutput		= & NewUserOutput;

//	dom.WriteXDMF("maz");
   	dom.Solve(/*tf*/4.2,/*dt*/t,/*dtOut*/0.01,"test",1000);

        return 0;
}
MECHSYS_CATCH
