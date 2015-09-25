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


using std::cout;
using std::endl;

void UserMazi(SPH::Domain & domi)
{
	double omega;
	omega = 1.16e-3;
	for (size_t a=0; a<domi.Particles.Size(); a++)
	{
		if (domi.Particles[a]->ID == 3)
		{
			domi.Particles[a]->translate(0.022);
			domi.Particles[a]->v  =-domi.Particles[a]->x(1)*omega,domi.Particles[a]->x(0)*omega,0.0;
			domi.Particles[a]->vb =-domi.Particles[a]->x(1)*omega,domi.Particles[a]->x(0)*omega,0.0;
		}
	}
}
//    	dom.GeneralBefore	= & UserMazi;
//    	dom.GeneralAfter	= & UserMazi;


int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Cs			= 3.0e-4;
//        dom.BC.Periodic[0]	= true;
        dom.Nproc		= 8;
    	dom.PresEq		= 0;
    	dom.VisEq		= 3;
    	dom.KernelType	= 4;
    	dom.NoSlip		= true;
//        dom.P0			= dom.Cs*dom.Cs*998.21*0.0005;
//    	dom.Shepard		= false;

        double xb,yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 5.0e-4;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

//    	dom.GeneralAfter	= & UserMazi;

        double maz;
        maz=(0.01*h/(dom.Cs));
        std::cout<<maz<<std::endl;

    	double R,Ri,Ro,no,mass;
    	Ri = 11.0e-3-3.5*dx;
    	Ro = 21.5e-3+3.5*dx;

    	for (size_t j=0;j<=((Ro-Ri)/dx);j++)
    	{
    		R = Ri+dx*j;
    		no = ceil(2*M_PI*R/dx);
    		for (size_t i=0; i<no; i++)
    		{
    			xb = R*cos(2*M_PI/no*i);
    			yb = R*sin(2*M_PI/no*i);
    			dom.AddSingleParticle(2,Vec3_t ( xb ,  yb , 0.0 ), 0.0 , rho , h , false);
    		}
     	}
		mass = M_PI*((Ro+dx/2.0)*(Ro+dx/2.0)-(Ri-dx/2.0)*(Ri-dx/2.0))*rho/dom.Particles.Size();

    	double omega;
    	omega = 1.16e-3;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->Mu = 1.002e-3;
    		dom.Particles[a]->MuRef = 1.002e-3;
    		dom.Particles[a]->Mass = mass;

    		xb=dom.Particles[a]->x(0);
    		yb=dom.Particles[a]->x(1);
    		if (sqrt(xb*xb+yb*yb)>21.5e-3)
    		{
    			R = sqrt(xb*xb+yb*yb);
    			dom.Particles[a]->ID		=3;
    			dom.Particles[a]->IsFree	=false;
    			dom.Particles[a]->v			=-yb/R*R*omega,xb/R*R*omega,0.0;
    			dom.Particles[a]->vb		=-yb/R*R*omega,xb/R*R*omega,0.0;
   		}
    		if (sqrt(xb*xb+yb*yb)<11.0e-3)
    		{
    			dom.Particles[a]->ID		=1;
    			dom.Particles[a]->IsFree	=false;
    		}
    	}

    	// No-Slip Particles
    	double m=2.0;
		R = 11.0e-3;
		no = ceil(2*M_PI*R/(m*dx));
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i);
			yb = R*sin(2*M_PI/no*i);
			dom.AddNSSingleParticle(1,Vec3_t ( xb ,  yb , 0.0 ), true);
		}
		R = 21.5e-3;
		no = ceil(2*M_PI*R/(m*dx));
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i);
			yb = R*sin(2*M_PI/no*i);
			dom.AddNSSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), true);
		}


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/2.0,"test06",1000);
        return 0;
}
MECHSYS_CATCH
