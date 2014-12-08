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
double a = 998.21*9.81*0.0547/(2.0*1.002e-3);
double b = 998.21*9.81*0.026/(2.0*1.002e-3);
double c = 998.21*9.81*0.006/(2.0*1.002e-3);

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(a*(2*0.001*position(1)-position(1)*position(1)),0.0,0.0);
	Den = 998.21*(1+9.81*(0.001-position(1))/(Cs*Cs));
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(0)<(5*h1))
	{
		Vel = Vec3_t(a*(2*0.001*position(1)-position(1)*position(1)),0.0,0.0);
		Den = 998.21*(1+9.81*(0.001-position(1))/(Cs*Cs));
	}
	else
	{
		Vel = Vec3_t(b*(2*0.001*position(1)-position(1)*position(1)),0.0,0.0);
		Den = 998.21*(1+9.81*(0.001-position(1))/(Cs*Cs));
	}
}
void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(c*(2*0.00209*position(1)-position(1)*position(1)),0.0,0.0);
	Den = 998.21*(1+9.81*(0.00209-position(1))/(Cs*Cs));
}

int main(int argc, char **argv) try
{
	SPH::Domain		dom;

	dom.Gravity		= 0.0,-9.81,0.0;
	dom.Dimension	= 2;
	dom.MU			= 1.002e-3;
	dom.Nproc		= 24;
	dom.PresEq		= 0;
	dom.VisEq		= 0;
	dom.KernelType	= 0;
	dom.Shepard		= true;
	dom.Cs			= 2.7;
	dom.TI			= 0.05;

	dom.BC.InOutFlow =3;
	dom.BC.allv = 0.0852,0.0,0.0;
	dom.BC.inv  = 0.1783,0.0,0.0;
	dom.BC.outv = 0.0852,0.0,0.0;
	dom.BC.inDensity = 998.21;
	dom.BC.outDensity = 998.21;
	dom.BC.allDensity = 998.21;

	dom.InCon = & UserInFlowCon;
	dom.AllCon = & UserAllFlowCon;
	dom.OutCon = & UserOutFlowCon;

	double yb,h,rho,H,U;
	double dx;

	rho	= 998.21;
	H	= 0.001;
	U	= 0.27;
	dx	= H/70.0;
	h	= dx*1.1;
	Cs	= dom.Cs;
	h1  = H;

	dom.InitialDist	= dx;
	dom.DomMax(1) = 2.5*H;

	double maz = (0.25*h/(dom.Cs+U));

	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -5.0*dx , 0.0 ), 15.0*H , H+5.0*dx  ,  0 , dx/2.0 ,rho, h, 0 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		yb=dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
		}
	}

	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/0.0005,"test06");
	return 0;
}
MECHSYS_CATCH
