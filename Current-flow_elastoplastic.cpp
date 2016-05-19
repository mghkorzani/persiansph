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

double H,U,RhoF,g,D,CsW,DampF,DampS,DampTime,dx;
int Check=0;

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=D)
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D) && position(1)>D)
		Vel = pow(((position(1)-D)/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D))
		Vel = 1.07*U , 0.0 , 0.0;

	Den = RhoF*pow((1+7.0*g*(H+D-position(1))/(CsW*CsW)),(1.0/7.0));
}

void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=D)
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D) && position(1)>D)
		Vel = pow(((position(1)-D)/(0.32*H)),(1.0/7.0))*0.92*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D))
		Vel = 1.07*0.92*U , 0.0 , 0.0;

	Den = RhoF*pow((1+7.0*g*(H+D-position(1))/(CsW*CsW)),(1.0/7.0));
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=D)
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D) && position(1)>D)
		Vel = pow(((position(1)-D)/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D))
		Vel = 1.07*U , 0.0 , 0.0;
}

void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<DampTime)
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
	  		if (domi.Particles[i]->ID == 2)
	        		domi.Particles[i]->a	= 0.0;
	}
	else
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
	  		if (domi.Particles[i]->ID == 2 && domi.Particles[i]->x(0)<0.4)
	        		domi.Particles[i]->a	= domi.Particles[i]->a*( 1.0-1.0/( 1+exp( 100.0*( domi.Particles[i]->x(0)-0.3 ) ) ) );
	  		if (domi.Particles[i]->ID == 2 && domi.Particles[i]->x(0)>1.0)
	        		domi.Particles[i]->a	= domi.Particles[i]->a*( 1.0/( 1+exp( 100.0*( domi.Particles[i]->x(0)-1.1) ) ) );
		}	
	}

//	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
//	for (size_t i=0; i<domi.Particles.Size(); i++)
//	{
//		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 1) domi.Particles[i]->a -= DampF * domi.Particles[i]->v;
//		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;
//	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
        SPH::Domain	dom;

        dom.Dimension	= 2;
	dom.SeepageType = 0;
        dom.Nproc	= 24;
    	dom.VisEq	= 3;
    	dom.KernelType	= 4;
    	dom.Scheme	= 0;
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;
	g		= norm(dom.Gravity);

    	double x,y,Muw,h,t,t1,t2,L;
    	dx	= 0.00625;
    	h	= dx*1.1;
	D	= 0.1;
	H	= 4.0*D;
	x	= 6.0*D;
	y	= 1.5*D+1.0*dx;
	L	= 13.0*D;

	U	= 0.4;
	RhoF	= 998.21;
	CsW	= 30.0;
	Muw	= 1.002e-3;
    	DampF	= 0.05*CsW/h;
        t1	= (0.25*h/(CsW));


    	DampTime		= 3.0;
    	dom.InitialDist 	= dx;
        dom.GeneralAfter	= & UserDamping;

	dom.BC.InOutFlow	= 3;
	dom.BC.inv		= U,0.0,0.0;
	dom.BC.outv		= U,0.0,0.0;
	dom.BC.allv		= U,0.0,0.0;
	dom.BC.inDensity	= RhoF;
	dom.BC.outDensity	= RhoF;
//	dom.BC.MassConservation = true;
	dom.InCon		= & UserInFlowCon;
	dom.OutCon		= & UserOutFlowCon;
	dom.AllCon		= & UserAllFlowCon;

    	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -4.0*dx , 0.0 ), L + dx/10.0 , 5.0*D + 8.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

    	double yb,xb,R,mass,no;;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		xb=dom.Particles[a]->x(0);
    		yb=dom.Particles[a]->x(1);

    		dom.Particles[a]->Cs		= CsW;
    		dom.Particles[a]->Alpha		= 0.05;
    		dom.Particles[a]->PresEq	= 1;
    		dom.Particles[a]->Mu		= Muw;
    		dom.Particles[a]->MuRef		= Muw;
    		dom.Particles[a]->Material	= 1;
    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H+D-yb)/(CsW*CsW)),(1.0/7.0));
    		dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H+D-yb)/(CsW*CsW)),(1.0/7.0));


    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    		if (yb>(H+D))
    		{
    			dom.Particles[a]->ID		= 6;
    			dom.Particles[a]->v		= 1.07*U , 0.0 , 0.0;
    			dom.Particles[a]->vb		= 1.07*U , 0.0 , 0.0;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    		if ((xb<=(5.0*dx) || xb>=(L-5.0*dx)) && yb<=(1.0*D))
    		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
       		if ((pow((xb-x),2)+pow((yb-y),2))<pow((D/2.0+dx/2.0),2.0))
       		{
        		dom.Particles[a]->ID		=4;
      		}
    	}

   	dom.DelParticles(4);
    	mass = RhoF*dx*dx;
	R = D/2.0-dx/2.0;
    	for (size_t j=0;j<5;j++)
    	{
    		if (j>0) R -= dx*0.85;
    		no = ceil(2*M_PI*R/dx);
    		for (size_t i=0; i<no; i++)
    		{
    			xb = x + R*cos(2*M_PI/no*i);
    			yb = y + R*sin(2*M_PI/no*i);
    			dom.AddSingleParticle(5,Vec3_t ( xb ,  yb , 0.0 ), mass , RhoF , h , true);
        		dom.Particles[dom.Particles.Size()-1]->Cs	= CsW;
        		dom.Particles[dom.Particles.Size()-1]->Alpha	= 0.05;
        		dom.Particles[dom.Particles.Size()-1]->PresEq	= 1;
        		dom.Particles[dom.Particles.Size()-1]->Mu	= 1.002e-3;
        		dom.Particles[dom.Particles.Size()-1]->MuRef	= 1.002e-3;
        		dom.Particles[dom.Particles.Size()-1]->Material	= 1;
        		dom.Particles[dom.Particles.Size()-1]->NoSlip	= true;
       		}
    	}

	double Nu,E,K,G,CsS,RhoS,c,Phi,Psi,d,n,k;

	Nu	= 0.25;
	E	= 10.0e6;
	K	= E/(3.0*(1.0-2.0*Nu));
	G	= E/(2.0*(1.0+Nu));
	n	= 0.5;
	RhoS	= 2650.0*(1.0-n)+n*RhoF;;
	CsS	= sqrt(K/RhoS);
//	CsS	= 80.0;
	c	= 0.0;
	Phi	= 180.0/M_PI * atan(0.389 / ( dx/2.0*g * (RhoS-RhoF) ) );
	Psi	= 0.0;
	d	= 0.00036;
	k	= n*n*n*d*d/(180.0*(1-n)*(1-n))*50.0;
  	DampS	= 0.02*sqrt(E/(RhoS*h*h));
        t2	= (0.25*h/(CsS));

        std::cout<<"CsS  = "<<CsS<<std::endl;
        std::cout<<"RhoS = "<<RhoS<<std::endl;
        std::cout<<"Phi  = "<<Phi<<std::endl;
        std::cout<<"K(Permeability) = "<<k<<std::endl;
        std::cout<<"K(Conductivity) = "<<k*RhoF*norm(dom.Gravity)/Muw<<std::endl;

	dom.AddBoxLength(2 ,Vec3_t ( 0.0 , -4.0*dx , 0.0 ), L + dx/10.0 , 1.0*D + 4.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, false );
//	dom.AddBoxLength(10,Vec3_t ( 12.0*D-15.0*dx , D , 0.0 ), 15.0*dx + dx/10.0 , 7.0*dx + dx/10.0     ,  0 , dx/2.0 ,RhoS, h, 1 , 0 , false, true );

    	mass = RhoS*dx*dx;
	R = D/2.0-dx/2.0;
    	for (size_t j=0;j<5;j++)
    	{
    		if (j>0) R -= dx*0.85;
    		no = ceil(2*M_PI*R/dx);
    		for (size_t i=0; i<no; i++)
    		{
    			xb = x + R*cos(2*M_PI/no*i);
    			yb = y + R*sin(2*M_PI/no*i);
    			dom.AddSingleParticle(10,Vec3_t ( xb ,  yb , 0.0 ), mass , RhoS , h , true);
//        		dom.Particles[dom.Particles.Size()-1]->NoSlip	= true;
       		}
    	}

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==2 || dom.Particles[a]->ID==10)
		{
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
//			dom.Particles[a]->TI		= 0.5;
//			dom.Particles[a]->TIn		= 2.55;
			dom.Particles[a]->d		= d;
			dom.Particles[a]->n		= n;
			dom.Particles[a]->k		= k;
			dom.Particles[a]->RhoF		= RhoF;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G		= G;
			dom.Particles[a]->K		= K;
			dom.Particles[a]->Fail		= 3;
			dom.Particles[a]->c		= c;
			dom.Particles[a]->phi		= Phi/180.0*M_PI;
			dom.Particles[a]->psi		= Psi/180.0*M_PI;

			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);

			if (dom.Particles[a]->ID == 10)
			{
				dom.Particles[a]->k = 100000.0;
			}
			if (yb<0.0)
			{
				dom.Particles[a]->ID	= 4;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= true;
			}
	    		if ((xb<=(5.0*dx) || xb>=(L-5.0*dx)) && yb<=(1.0*D))
    			{
				dom.Particles[a]->ID	= 4;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= true;
			}

		}
	}

        t	= std::min(t1,t2);
        std::cout<<"t1 = "<<t1<<std::endl;
        std::cout<<"t2 = "<<t2<<std::endl;
        std::cout<<"t  = "<<t<<std::endl;


   	dom.Solve(/*tf*/700.0,/*dt*/t,/*dtOut*/0.1,"test",100000);
        return 0;
}
MECHSYS_CATCH
