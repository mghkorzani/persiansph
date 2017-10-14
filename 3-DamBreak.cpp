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

int main(int argc, char **argv) try
{
  SPH::Domain		dom;

  dom.Dimension   = 2;
  dom.Nproc       = 4;
  dom.Scheme			= 0;
	dom.Viscosity_Eq_Set(Morris);
	dom.Kernel_Set(Qubic_Spline);
	dom.Gradient_Approach_Set(Squared_density);

  double xb,yb,dx,Cs,h,Rho,H,L,TL,TH,t,res,Mu,g;

  H   = 0.6;
  L   = 2.0*H;
  TH  = 1.0;
  TL  = 5.366*H;
  res = 40.0;

  g   = 9.81;
  Rho = 998.21;
  Mu	= 1.002e-3;
  dx  = H/res;
  h   = dx*1.2;
  Cs	= 20.0 * sqrt(g*H);
  t   = (0.2*h/Cs);

  dom.InitialDist 	= dx;
  dom.Gravity	      = 0.0, -g ,0.0 ;

	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , 0.0 ), 7.0*dx + TL + dx/10.0 , 3.0*dx + TH + dx/10.0,  0 , dx/2.0 ,Rho, h, 1 , 0 , false, false );

  for (size_t a=0; a<dom.Particles.Size(); a++)
  {
    xb=dom.Particles[a]->x(0);
    yb=dom.Particles[a]->x(1);

    dom.Particles[a]->Cs      = Cs;
		dom.Particles[a]->PresEq	= 1;
		dom.Particles[a]->Mu			= Mu;
		dom.Particles[a]->MuRef		= Mu;
		dom.Particles[a]->Material= 1;
    dom.Particles[a]->Shepard = false;


    if (xb<0.0)
    {
			dom.Particles[a]->ID			= 2;
			dom.Particles[a]->IsFree	= false;
			dom.Particles[a]->NoSlip	= true;
    }

    if (yb<0.0)
    {
			dom.Particles[a]->ID			= 2;
			dom.Particles[a]->IsFree	= false;
			dom.Particles[a]->NoSlip	= true;
    }

    if (xb>TL)
    {
			dom.Particles[a]->ID			= 2;
			dom.Particles[a]->IsFree	= false;
			dom.Particles[a]->NoSlip	= true;
    }

    if (yb>H && dom.Particles[a]->ID==1)
      dom.Particles[a]->ID=11;

    if (xb>L && dom.Particles[a]->ID==1)
      dom.Particles[a]->ID=11;

    if (dom.Particles[a]->ID==1)
    {
      dom.Particles[a]->Density	  = Rho*pow((1+7.0*g*(H-yb)/(Cs*Cs)),(1.0/7.0));
      dom.Particles[a]->Densityb	= Rho*pow((1+7.0*g*(H-yb)/(Cs*Cs)),(1.0/7.0));
    }

  }

  dom.DelParticles(11);

  dom.Solve(/*tf*/50.0,/*dt*/t,/*dtOut*/0.05,"test06",999);
  return 0;
}
MECHSYS_CATCH
