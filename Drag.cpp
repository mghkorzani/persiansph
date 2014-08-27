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

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.MU			= 1.002e-3;
	dom.P0			= 0.0003;
	dom.PresEq		= 0;
	dom.VisEq		= 0;
	dom.KernelType	= 2;
	dom.Nproc		= 24;

	dom.TI			= 0.05;

	double xb,yb,h,rho,mass;
	double dx,R,Rc,Re;
	size_t no;

	rho = 998.21;
	h = 0.002*1.2;
	dx = 0.002;
	Rc = 0.02;
	mass = (sqrt(3.0)*dx*dx/4.0)*rho;
	Re = 4.0;

	dom.ConstVelPeriodic= Re*dom.MU/(rho*2.0*Rc);
	dom.Cs				= dom.ConstVelPeriodic*20.0;
	dom.InitialDist = dx;

	dom.AddRandomBox(3 ,Vec3_t ( -100.0*0.002 , -100.0*0.002 , 0.0 ), 200.0*0.002 ,200.0*0.002  ,  0 , 0.001 ,rho, h);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb)<((Rc+3.0/2.0*sqrt(3.0)*dx)*(Rc+3.0/2.0*sqrt(3.0)*dx)))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	R = Rc+sqrt(3.0)*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , false);
	}

	R = Rc+sqrt(3.0)/2.0*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , false);
	}

	R = Rc;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
	}

	R = Rc-sqrt(3.0)/2.0*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
	}

	R = R-sqrt(3.0)/2*dx;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
	}

	dom.Solve(/*tf*/10000.0,/*dt*/(0.2*h/(dom.Cs+dom.ConstVelPeriodic)),/*dtOut*/(2.0*h/(dom.Cs+dom.ConstVelPeriodic)),"test06");
	return 0;
}
MECHSYS_CATCH
