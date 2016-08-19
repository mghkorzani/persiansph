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
	double g,dx,RhoF,CsW,DampTime,DampF,DampS,L,H,V;
	int	check=0;

void UserInFlowCon(SPH::Domain & domi)
{
	if (domi.Time>=DampTime)
	{
		int temp;
		for (int q1=0;  q1<2                ; q1++)
		for (int q2=0;  q2<domi.CellNo[1]   ; q2++)
		for (int q3=0;  q3<domi.CellNo[2]   ; q3++)
			if (domi.HOC[q1][q2][q3]!=-1)
			{
				temp = domi.HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (domi.Particles[temp]->IsFree)
					{
						domi.Particles[temp]->dDensity	= 0.0;
						domi.Particles[temp]->Density	= RhoF*pow((1+7.0*(g*(H-domi.Particles[temp]->x(1))+0.5*V*V)/(CsW*CsW)),(1.0/7.0));
						domi.Particles[temp]->Densityb	= RhoF*pow((1+7.0*(g*(H-domi.Particles[temp]->x(1))+0.5*V*V)/(CsW*CsW)),(1.0/7.0));
		    				domi.Particles[temp]->Pressure	= SPH::EOS(domi.Particles[temp]->PresEq, domi.Particles[temp]->Cs, domi.Particles[temp]->P0,
													domi.Particles[temp]->Density,domi.Particles[temp]->RefDensity);
						domi.Particles[temp]->Mu	= domi.Particles[temp]->MuRef; 	
						domi.Particles[temp]->v  	= 0.0 , V , 0.0;
						domi.Particles[temp]->vb 	= 0.0 , V , 0.0;
					}
					temp = domi.Particles[temp]->LL;
				}
			}
	}
}


void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<DampTime)
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
		if (check==0)
		{
			#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
			for (size_t i=0; i<domi.Particles.Size(); i++)
			{
				if (domi.Particles[i]->IsFree)
				{ 
					domi.Particles[i]->v = 0.0;;
					domi.Particles[i]->vb = 0.0;;
				}
			}
			domi.SWIType	= 0;
		}
		check = 1;

		int temp;
		for (int q1=0;  q1<<domi.CellNo[0]  ; q1++)
		for (int q2=0;  q2<2                ; q2++)
		for (int q3=0;  q3<domi.CellNo[2]   ; q3++)
			if (domi.HOC[q1][q2][q3]!=-1)
			{
				temp = domi.HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (domi.Particles[temp]->IsFree)
					{
						domi.Particles[temp]->dDensity  = 0.0;
						domi.Particles[temp]->a		= 0.0 , 0.0 , 0.0;
					}
					temp = domi.Particles[temp]->LL;
				}
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
    	dom.KernelType	= 4;
	dom.SWIType	= 3;
    	dom.Scheme	= 0;
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;
	g		= norm(dom.Gravity);

        dom.BC.Periodic[1]	= true;
        dom.BC.Periodic[0]	= true;

        dom.GeneralBefore	= & UserInFlowCon;
        dom.GeneralAfter	= & UserDamping;


    	double h,t,t1,t2,Muw,Q;
    	dx	= 0.002;
    	h	= dx*1.2;
	dom.InitialDist	= dx;

//	L	= 0.15;
	L	= 0.05;
	H	= 0.25;
	Q	= 0.126;
	V	= (Q/1000.0)/(0.15*0.15);
	RhoF	= 1000.0;
	CsW	= 10.0*sqrt(2.0*g*H)*2.0;
	Muw	= 1.0e-3;
        t1	= (0.25*h/(CsW));



//    	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-4.0*dx , 0.0 , 0.0 ), L + 8.0*dx + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );
    	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , 0.0 , 0.0 ), L + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

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
    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H-yb)/(CsW*CsW)),(1.0/7.0));
		dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H-yb)/(CsW*CsW)),(1.0/7.0));
/*
    		if (xb<-L/2.0 || xb>L/2.0)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
*/
   	}

	double Nu,E,K,G,CsS,RhoS,c,Phi,Psi,n1,n2,n3,d1,d2,d3,h1,h2,h3,h0;

	Nu	= 0.3;
	E	= 25.0e6;
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
        t2	= (0.2*h/(CsS));
	h0	= 0.05;
	h1	= 0.025;
	h2	= 0.025;
	h3	= 0.11;

        std::cout<<"t2 = "<<t2<<std::endl;

        std::cout<<"CsS  = "<<CsS<<std::endl;
        std::cout<<"RhoS = "<<RhoS<<std::endl;
        std::cout<<"Phi  = "<<Phi<<std::endl;
        std::cout<<"C    = "<<c<<std::endl;

//	dom.AddBoxLength(5 ,Vec3_t ( -L/2.0-4.0*dx , h0 , 0.0 ), L + 8.0*dx + dx/10.0 , H-h0 + dx/10.0 ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, false );
	dom.AddBoxLength(5 ,Vec3_t ( -L/2.0 , h0 , 0.0 ), L + dx/10.0 , H-h0 + dx/10.0 ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==5)
		{
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Alpha		= 0.2;
			dom.Particles[a]->Beta		= 0.2;
			if (c>0.0)
			{
				dom.Particles[a]->TI	= 0.5;
				dom.Particles[a]->TIn	= 2.55;
			}
			dom.Particles[a]->TIInitDist	= dx;
			dom.Particles[a]->VarPorosity	= true;
			dom.Particles[a]->SeepageType	= 2;	
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
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
				dom.Particles[a]->d		= d1*1.0e20;
	    		}
*/
	    		if (yb>(h0+h1+h2+h3) && dom.Particles[a]->ID==5)
	    			dom.Particles[a]->ID		= 10;

	    		if (yb>(h0) && yb<=(h0+h1) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 5;
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
				dom.Particles[a]->d		= d1;
				dom.Particles[a]->n0		= n1;
	    		}
	    		if (yb>(h0+h1) && yb<=(h0+h1+h2) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 6;
	    			dom.Particles[a]->IsFree	= false;
	    			dom.Particles[a]->NoSlip	= false;
				dom.Particles[a]->d		= d2;
				dom.Particles[a]->n0		= n2;
	    		}
	    		if (yb>(h0+h1+h2) && yb<=(h0+h1+h2+h3) && dom.Particles[a]->ID==5)
	    		{
	    			dom.Particles[a]->ID		= 7;
				dom.Particles[a]->d		= d3;
				dom.Particles[a]->n0		= n3;
	    		}
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
	dom.OutputName[2]	= "Permeability";
        dom.UserOutput		= & NewUserOutput;

   	dom.Solve(/*tf*/700.0,/*dt*/t,/*dtOut*/0.01,"test",100000);
        return 0;
}
MECHSYS_CATCH
