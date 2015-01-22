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

double Cs,h1;

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(1.5,0.0,0.0);
	Den = 998.21*(1+9.81*(0.036-position(1))/(Cs*Cs));
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(0)<1.6)
	{
		Vel = Vec3_t(1.5,0.0,0.0);
		Den = 998.21*(1+9.81*(0.036-position(1))/(Cs*Cs));
	}
	else
	{
		Vel = Vec3_t(0.46,0.0,0.0);
		Den = 998.21*(1+9.81*(0.036-position(1))/(Cs*Cs));
	}
}
//void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
//{
//	Vel = Vec3_t(0.46,0.0,0.0);
//	Den = 998.21*(1+9.81*(0.036-position(1))/(Cs*Cs));
//}

int main(int argc, char **argv) try
{
	double Tail;
	Tail = atof(argv[1]);

	SPH::Domain		dom;

	dom.Gravity		= 0.0,-9.81,0.0;
	dom.Dimension	= 2;
	dom.MU			= 1.002e-3;
	dom.Nproc		= 20;
	dom.PresEq		= 0;
	dom.VisEq		= 1;
	dom.KernelType	= 4;
	dom.Shepard		= true;
	dom.TI			= 0.025;

	dom.BC.InOutFlow =3;
	dom.BC.allv = 1.5,0.0,0.0;
	dom.BC.inv  = 1.5,0.0,0.0;
//	dom.BC.outv = 0.46,0.0,0.0;
	dom.BC.inDensity = 998.21;
//	dom.BC.outDensity = 998.21;
	dom.BC.allDensity = 998.21;

	dom.InCon  = & UserInFlowCon;
	dom.AllCon = & UserAllFlowCon;
//	dom.OutCon = & UserOutFlowCon;

	double xb,yb,h,rho,H,U;
	double dx;

	rho	= 998.21;
	H	= 0.036;
	U	= 2.0;
	dom.Cs = 10.0*U;
	dx	= H/9.0;
	h	= dx*1.1;
	Cs	= dom.Cs;
	h1  = H;

	dom.InitialDist	= dx;
//	dom.DomMax(1) = 5.0*H;

	double maz = (0.15*h/(dom.Cs+U));

	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -5.0*dx , 0.0 ), 3.2 , 5.0*H+5.0*dx  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->Density = 998.21*(1+9.81*(0.1-dom.Particles[a]->x(1))/(Cs*Cs));
			dom.Particles[a]->Densityb = 998.21*(1+9.81*(0.1-dom.Particles[a]->x(1))/(Cs*Cs));
			dom.Particles[a]->IsFree=false;
		}
		if (xb>(3.2-5.0*dx) && yb<Tail)
		{
			dom.Particles[a]->ID=3;
			dom.Particles[a]->IsFree=false;
		}
		if (xb<=H && yb>H)
		{
			dom.Particles[a]->ID=3;
			dom.Particles[a]->IsFree=false;
		}
		if (xb>H && yb>H)
		{
			dom.Particles[a]->ID=4;
		}
	}
	dom.DelParticles(4);

	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/0.1,"test06",2);
	return 0;
}
MECHSYS_CATCH
