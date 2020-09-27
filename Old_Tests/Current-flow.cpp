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

double H,U,RhoF,RhoS,g,D,CsW,CsS,DampTime,dx,T0,L,Muw;
double Z0,U0,Z,Us;

void UserInFlowCon(SPH::Domain & domi)
{
	Us = U0/7.0*pow((3.6e-4/(H+D-Z0)),(1.0/7.0));
	Z = 2.5*3.6e-4/30.0*(1.0-exp(-Us*2.5*3.6e-4/(27.0*Muw/RhoF))) + (Muw/RhoF)/(9.0*Us);

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
					if (domi.Particles[temp]->ID==1 && domi.Particles[temp]->x(1)<= (H+D))
					{
						
						domi.Particles[temp]->dDensity = 0.0;
						domi.Particles[temp]->Density  = RhoF*pow((1+7.0*g*(H+D-domi.Particles[temp]->x(1))/(CsW*CsW)),(1.0/7.0));
						domi.Particles[temp]->Densityb = RhoF*pow((1+7.0*g*(H+D-domi.Particles[temp]->x(1))/(CsW*CsW)),(1.0/7.0));
	    					domi.Particles[temp]->Pressure = SPH::EOS(domi.Particles[temp]->PresEq, domi.Particles[temp]->Cs, domi.Particles[temp]->P0,
													domi.Particles[temp]->Density,domi.Particles[temp]->RefDensity);	
					}
					if (domi.Particles[temp]->ID==2 && domi.Particles[temp]->x(1)<= D)
					{
						domi.Particles[temp]->dDensity = 0.0;
						domi.Particles[temp]->Density  = RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*H+RhoS*(D-domi.Particles[temp]->x(1)))),(1.0/7.0));
						domi.Particles[temp]->Densityb = RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*H+RhoS*(D-domi.Particles[temp]->x(1)))),(1.0/7.0));
	    					domi.Particles[temp]->Pressure = SPH::EOS(domi.Particles[temp]->PresEq, domi.Particles[temp]->Cs, domi.Particles[temp]->P0,
													domi.Particles[temp]->Density,domi.Particles[temp]->RefDensity);	
					}
					if (domi.Particles[temp]->x(1)> Z0)
					{
						domi.Particles[temp]->v  = Us/0.4*log((domi.Particles[temp]->x(1)-Z0)/Z) , 0.0 , 0.0;
						domi.Particles[temp]->vb = Us/0.4*log((domi.Particles[temp]->x(1)-Z0)/Z) , 0.0 , 0.0;
					}
					if (domi.Particles[temp]->x(1)<= Z0)
					{
						domi.Particles[temp]->v(1)  = 0.0;
						domi.Particles[temp]->vb(1) = 0.0;
						if (domi.Particles[temp]->v(0)<0 )
						{
							domi.Particles[temp]->v(0) = 0.0;
							domi.Particles[temp]->vb(0) = 0.0;
						}
					}
				}
				temp = domi.Particles[temp]->LL;
			}
		}
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Us = U0/7.0*pow((3.6e-4/(H+D-Z0)),(1.0/7.0));
	Z = 2.5*3.6e-4/30.0*(1.0-exp(-Us*2.5*3.6e-4/(27.0*Muw/RhoF))) + (Muw/RhoF)/(9.0*Us);

	if (position(1)>Z0)
		Vel = Us/0.4*log((position(1)-Z0)/Z) , 0.0 , 0.0;
}

void UserDamping(SPH::Domain & domi)
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
					if (domi.Particles[temp]->ID==1  && domi.Particles[temp]->x(1)<= (H+D))
						domi.Particles[temp]->dDensity  = 0.0;
					if (domi.Particles[temp]->ID==2 && domi.Particles[temp]->x(1)<= D)
						domi.Particles[temp]->dDensity  = 0.0;
					if (domi.Particles[temp]->x(1)>Z0)
					{
						domi.Particles[temp]->a = 0.0 , 0.0 , 0.0;
					}
					if (domi.Particles[temp]->x(1)<=Z0 && domi.Particles[temp]->ID==1)
					{
						domi.Particles[temp]->a(1) = 0.0;
					}
				}
				temp = domi.Particles[temp]->LL;
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
	double IU;
	infile >> IT0;	infile.ignore(200,'\n');
	infile >> Im;	infile.ignore(200,'\n');
	infile >> IU;	infile.ignore(200,'\n');
        SPH::Domain	dom;

        dom.Dimension	= 2;
        dom.Nproc	= 24;
    	dom.VisEq	= 3;
    	dom.KernelType	= 4;
    	dom.Scheme	= 0;
    	dom.Gravity	= 0.0 , -9.81 , 0.0 ;


    	double x,y,Mus,h,t,t1,t2,t3,m,n;
	// Geometrical parameters
    	dx	= 0.01;
    	dom.InitialDist = dx;
    	h	= dx*1.1;
	D	= 0.1;
	g	= norm(dom.Gravity);
	H	= 4.0*D;
	L	= 10.0*D;
	U	= IU;


	x	= 4.0*D;
	y	= 1.5*D + 1.1*dx;
	Z0	= D + D*3.0/4.0;
	U0	= 0.342/0.4 * U;

	// Water parameters
	RhoF	= 998.21;
	CsW	= 25.0;
	Muw	= 1.002e-3;

	// Sediment parameters
	CsS	= 35.0;
	T0	= IT0;
	m	= Im;
	Mus	= 0.05;
	n	= 0.5;
	RhoS	= 2650.0*(1.0-n)+n*RhoF;
	RhoS	= 1300.0;

        std::cout<<"Uavg = "<<U<<std::endl;
        std::cout<<"Ty   = "<<T0<<std::endl;
        std::cout<<"m    = "<<m<<std::endl;
        std::cout<<"MuB  = "<<Mus<<std::endl;
        std::cout<<"CsS  = "<<CsS<<std::endl;
        std::cout<<"RhoS = "<<RhoS<<std::endl;

	// Time step selection
        t1	= (0.25*h/(CsW));
        t2	= (0.25*h/(CsS));
        t3	= 0.5*0.125*h*h*(RhoS-RhoF)/(Mus+m*T0);
        t	= std::min(t1,t2);
        t	= std::min(t,t3);

        std::cout<<"tw   = "<<t1<<std::endl;
        std::cout<<"ts   = "<<t2<<std::endl;
        std::cout<<"tvis = "<<t3<<std::endl;
        std::cout<<"tmin = "<<t<<std::endl;

	dom.BC.allv		= U,0.0,0.0;
	dom.AllCon		= & UserAllFlowCon;
        dom.GeneralBefore	= & UserInFlowCon;
        dom.GeneralAfter	= & UserDamping;
        dom.BC.Periodic[0]	= true;
	dom.DomMax(1)		= H+D+D/2.0;

//	Us = U0/7.0*pow((3.6e-4/(H+D-Z0)),(1.0/7.0));
//	Z = 2.5*3.6e-4/30.0*(1.0-exp(-Us*2.5*3.6e-4/(27.0*Muw/RhoF))) + (Muw/RhoF)/(9.0*Us);

    	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -4.0*dx , 0.0 ), L + dx/10.0 , 5.0*D + 4.0*dx + dx/10.0 ,  0 , dx/2.0 ,RhoF, h, 1 , 0 , false, false );

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
    		dom.Particles[a]->Density	= RhoF*pow((1+7.0*g*(H+D-yb)/(CsW*CsW)),(1.0/7.0));
    		dom.Particles[a]->Densityb	= RhoF*pow((1+7.0*g*(H+D-yb)/(CsW*CsW)),(1.0/7.0));


    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
/*    		if (yb>(H+D))
    		{
    			dom.Particles[a]->ID		= 3;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
			dom.Particles[a]->v  = Us/0.4*log((H+D-Z0)/Z) , 0.0 , 0.0;
			dom.Particles[a]->vb = Us/0.4*log((H+D-Z0)/Z) , 0.0 , 0.0;

    		}
*/     		if (yb>=0.0 && yb<=(1.0*D) && dom.Particles[a]->ID == 1)
    		{
    			dom.Particles[a]->ID		= 2;
        		dom.Particles[a]->Cs		= CsS;
        		dom.Particles[a]->Mu		= Mus;
        		dom.Particles[a]->MuRef		= Mus;
        		dom.Particles[a]->m 		= m;
        		dom.Particles[a]->T0 		= T0;
        		dom.Particles[a]->Density	= RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*H+RhoS*(D-yb))),(1.0/7.0));
        		dom.Particles[a]->Densityb	= RhoS*pow((1+7.0*g/(CsS*CsS*RhoS)*(RhoF*H+RhoS*(D-yb))),(1.0/7.0));
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
	R = D/2.0-dx/2.0;
    	for (size_t j=0;j<5;j++)
    	{
    		if (j>0) R -= dx*1.0;
    		no = ceil(2*M_PI*R/dx);
    		for (size_t i=0; i<no; i++)
    		{
    			xb = x + R*cos(2*M_PI/no*i);
    			yb = y + R*sin(2*M_PI/no*i);
    			dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), mass , RhoF , h , true);
        		dom.Particles[dom.Particles.Size()-1]->Cs	= CsW;
        		dom.Particles[dom.Particles.Size()-1]->Alpha	= 0.05;
        		dom.Particles[dom.Particles.Size()-1]->PresEq	= 1;
        		dom.Particles[dom.Particles.Size()-1]->Mu	= Muw;
        		dom.Particles[dom.Particles.Size()-1]->MuRef	= Muw;
        		dom.Particles[dom.Particles.Size()-1]->Material	= 1;
        		dom.Particles[dom.Particles.Size()-1]->NoSlip	= true;
       		}
    	}


   	dom.Solve(/*tf*/700.0,/*dt*/t,/*dtOut*/0.1,"Bin_Sed",10000);
        return 0;
}
MECHSYS_CATCH
