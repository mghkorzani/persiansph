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

#include <Domain.h>


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.NoSlip		= true;
	dom.PeriodicX	= true;
	dom.PeriodicY	= true;
	dom.ConstVelPeriodic= 4.461319*10.0;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= dom.ConstVelPeriodic*20.0;
	dom.MU			= 1.002e-3;
	dom.P0			= 63;
	dom.PresEq		= 0;
	dom.VisEq		= 0;
	dom.KernelType	= 2;
	dom.Nproc		= 24;

	dom.TI			= 0.05;
	dom.InitialDist = 0.01;

	double xb,yb;

	dom.AddRandomBox(3 ,Vec3_t ( -125*0.01 , -125*0.01 , 0.0 ), 250.0*0.01 ,250*0.01  ,  0 , 0.005 ,9.9821e-4, 0.012);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb<0.037))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	double dx,R,Rc;
	size_t no;
	dx = 0.01;

	Rc = 0.1125*1.5;


	R = Rc+sqrt(3.0)*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
	}

	R = Rc+sqrt(3.0)/2.0*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
	}

	R = Rc;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = Rc-sqrt(3.0)/2.0*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	dom.Solve(/*tf*/0.2,/*dt*/(0.2*0.012/(dom.Cs+dom.ConstVelPeriodic)),/*dtOut*/5.0e-5,"test06");
	return 0;
}
MECHSYS_CATCH
