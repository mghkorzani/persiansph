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

double H,U,U0,RhoF,RhoS,g,D,CsW,CsS,Damp,DampTime,dx,T0;
int Check=0;

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=(D+dx))
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D+dx) && position(1)>(D+dx))
		Vel = pow(((position(1)-(D+dx))/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D+dx))
		Vel = 1.07*U , 0.0 , 0.0;

	if (position(1)<=D)
		Den = RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*(H+dx)+RhoS*(D-position(1)))),(1.0/7.0));
	else
		Den = RhoF*pow((1+7.0*g*(H+D+dx-position(1))/(CsW*CsW)),(1.0/7.0));
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=(D+dx))
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D+dx) && position(1)>(D+dx))
		Vel = pow(((position(1)-(D+dx))/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D+dx))
		Vel = 1.07*U , 0.0 , 0.0;

	if (position(1)<=D)
		Den = RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*(H+dx)+RhoS*(D-position(1)))),(1.0/7.0));
	else
		Den = RhoF*pow((1+7.0*g*(H+D+dx-position(1))/(CsW*CsW)),(1.0/7.0));
}

void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(1)<=(D+dx+4.0*dx))
		Vel = 0.0 , 0.0 , 0.0;
	if (position(1)<(0.514*H+D+dx) && position(1)>(D+dx+4.0*dx))
		Vel = pow(((position(1)-(D+dx))/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
	if (position(1)>=(0.514*H+D+dx))
		Vel = 1.07*U , 0.0 , 0.0;

	if (position(1)<=D)
		Den = RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*(H+dx)+RhoS*(D-position(1)))),(1.0/7.0));
	else
		Den = RhoF*pow((1+7.0*g*(H+D+dx-position(1))/(CsW*CsW)),(1.0/7.0));
}

void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<DampTime)
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->IsFree) domi.Particles[i]->a -= Damp * domi.Particles[i]->v;

	}

	if (domi.Time>DampTime && Check==0)
	{
		domi.BC.inv	= U,0.0,0.0;
		domi.BC.outv	= U,0.0,0.0;
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
	  		if (domi.Particles[i]->ID == 2)
	        		domi.Particles[i]->T0	= T0;
			if (domi.Particles[i]->IsFree)
			{
				domi.Particles[i]->FirstStep = true;
				if (domi.Particles[i]->x(1)<=(D+dx))
				{
					domi.Particles[i]->v = 0.0 , 0.0 , 0.0;
					domi.Particles[i]->vb = 0.0 , 0.0 , 0.0;
				}
				if (domi.Particles[i]->x(1)<(0.514*H+D+dx) && domi.Particles[i]->x(1)>(D+dx))
				{
					domi.Particles[i]->v = pow(((domi.Particles[i]->x(1)-(D+dx))/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
					domi.Particles[i]->vb = pow(((domi.Particles[i]->x(1)-(D+dx))/(0.32*H)),(1.0/7.0))*U , 0.0 , 0.0;
				}
				if (domi.Particles[i]->x(1)>=(0.514*H+D+dx))
				{
					domi.Particles[i]->v = 1.07*U , 0.0 , 0.0;
					domi.Particles[i]->vb = 1.07*U , 0.0 , 0.0;
				}
			}
		}
		Check = 1;
	}
	
	if (domi.Time>DampTime && Check==1)
	{
		if (U<=U0)
			U = ((U0-0.01)/2.5)*(domi.Time-DampTime)+0.01;
		else
		{
			U = U0;
			Check = 2;
		}
	}
}


using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char **argv) try
{
	if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
	String filekey  (argv[1]);
	String filename (filekey+".inp");
	ifstream infile(filename.CStr());
	double IT0;
	double Im;
	infile >> IT0;	infile.ignore(200,'\n');
	infile >> Im;	infile.ignore(200,'\n');
        SPH::Domain	dom;

        dom.Dimension	= 2;
        dom.Nproc	= 10;
    	dom.VisEq	= 0;
    	dom.KernelType	= 4;
    	dom.Scheme	= 0;
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;
//    	dom.XSPH	= 0.5;

    	double x,y,Muw,Mus,h,t,t1,t2,t3,m,por;
    	dx	= 0.004;
    	h	= dx*1.1;
	D	= 0.1;
	g	= norm(dom.Gravity);
	H	= 4.0*D;
	U	= 0.01;
	U0	= 0.4;
	RhoF	= 998.21;
	CsW	= 30.0;
	x	= 6.0*D;
	y	= 1.5*D + dx;
	Muw	= 1.002e-3;
	Mus	= 0.05;
	CsS	= 35.0;
	T0	= IT0;
	m	= 5.0;

        std::cout<<"T0 = "<<T0<<std::endl;
        std::cout<<"m  = "<<m<<std::endl;

	por	= 0.4;
	RhoS	= 2650.0*(1.0-por)+por*RhoF;
        std::cout<<"RhoS = " <<RhoS<<std::endl;

        t1	= (0.25*h/(CsW));
        t2	= (0.25*h/(CsS));
        t3	= 0.125*h*h*(RhoS-RhoF)/(Mus+m*T0);
        t	= std::min(t1,t2);
        t	= std::min(t,t3);

        std::cout<<"t1 = "<<t1<<std::endl;
        std::cout<<"t2 = "<<t2<<std::endl;
        std::cout<<"t3 = "<<t3<<std::endl;
        std::cout<<"t  = "<<t<<std::endl;

    	Damp	 		= 0.05*CsW/h;
    	DampTime		= 0.5;
    	dom.InitialDist 	= dx;

	dom.BC.InOutFlow	= 3;
//	dom.BC.inv		= 0.000001,0.0,0.0;
//	dom.BC.out		= 0.000001,0.0,0.0;
	dom.BC.inDensity	= RhoF;
	dom.BC.outDensity	= RhoF;
//	dom.BC.MassConservation = true;

	dom.InCon		= & UserInFlowCon;
	dom.OutCon		= & UserOutFlowCon;
        dom.GeneralBefore	= & UserDamping;


    	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -4.0*dx , 0.0 ), 12.0*D + dx/10.0 , 5.0*D + 5.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

    	double yb,xb;

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
//    		dom.Particles[a]->Shepard	= true;
    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H+D+dx-yb)/(CsW*CsW)),(1.0/7.0));
    		dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H+D+dx-yb)/(CsW*CsW)),(1.0/7.0));


    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}

    		if ((xb<=(5.0*dx) || xb>=(12.0*D-5.0*dx)) && yb<=(1.0*D))
    		{
    			dom.Particles[a]->ID		= 5;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}

    		if (yb>=0.0 && yb<=(1.0*D) && dom.Particles[a]->ID == 1)
    		{
    			dom.Particles[a]->ID		= 2;
        		dom.Particles[a]->Cs		= CsS;
        		dom.Particles[a]->Mu		= Mus;
        		dom.Particles[a]->MuRef		= Mus;
        		dom.Particles[a]->m 		= m;
//        		dom.Particles[a]->T0 		= T0;
        		dom.Particles[a]->Density	= RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*(H+dx)+RhoS*(D-yb))),(1.0/7.0));
        		dom.Particles[a]->Densityb	= RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*(H+dx)+RhoS*(D-yb))),(1.0/7.0));
        		dom.Particles[a]->RefDensity	= RhoS;
        		dom.Particles[a]->Mass		= dom.Particles[a]->Mass/RhoF*RhoS;
    		}

       		if ((pow((xb-x),2)+pow((yb-y),2))<pow((D/2.0+dx/2.0),2.0))
       		{
        		dom.Particles[a]->ID=4;
       		}
    	}

    	dom.DelParticles(4);

    	double R,mass,no;
    	mass = RhoF*dx*dx;

    	for (size_t j=0;j<4;j++)
    	{
    		R = D/2.0-dx*j;
    		no = ceil(2*M_PI*R/dx);
    		for (size_t i=0; i<no; i++)
    		{
    			xb = x + R*cos(2*M_PI/no*i);
    			yb = y + R*sin(2*M_PI/no*i);
    			dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , RhoF , h , true);
        		dom.Particles[dom.Particles.Size()-1]->Cs	= CsW;
        		dom.Particles[dom.Particles.Size()-1]->Alpha	= 0.05;
        		dom.Particles[dom.Particles.Size()-1]->PresEq	= 1;
        		dom.Particles[dom.Particles.Size()-1]->Mu	= 1.002e-3;
        		dom.Particles[dom.Particles.Size()-1]->MuRef	= 1.002e-3;
        		dom.Particles[dom.Particles.Size()-1]->Material	= 1;
        		dom.Particles[dom.Particles.Size()-1]->NoSlip	= true;
       		}
    	}


   	dom.Solve(/*tf*/700.0,/*dt*/t,/*dtOut*/0.1,"Bin_Sed",10000);
        return 0;
}
MECHSYS_CATCH
