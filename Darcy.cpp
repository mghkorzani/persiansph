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
#include "./Source/Interaction.h"


using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.KernelType	= 0;
	dom.Nproc		= 24;

	dom.Scheme		= 0;
	dom.Gravity		= 0.0,-9.81,0.0;
	dom.BC.Periodic[0] = true;
	dom.BC.Periodic[1] = true;

	double h,rhoF,rhoS,CsF;
	double dx,HF,HS,L,T,k,n;

	HF		= 3.0;
	HS		= 3.0;
	L		= 0.6;
	rhoF	= 998.21;
	rhoS	= 2038.7;
	dx		= 0.05;
	h		= dx*1.3;
	CsF		= 10.0*sqrt(1.0*9.81*HF);
	k		= 5.0e-2;
	n		= 1.0;
	T		= (0.1*h/CsF);
	dom.InitialDist	= dx;
    std::cout<<"T  = "<<T<<std::endl;

	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , 0.0 , 0.0 ), L + dx/10.0 , HF + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
 		dom.Particles[a]->PresEq	= 0;
 		dom.Particles[a]->Alpha		= 0.03;
 		dom.Particles[a]->Mu		= 1.002e-3;
		dom.Particles[a]->MuRef		= 1.002e-3;
		dom.Particles[a]->Material	= 1;
		dom.Particles[a]->Cs		= CsF;
	}


	dom.AddBoxLength(3 ,Vec3_t ( -L/2.0 , 0.0 , 0.0 ), L + dx/10.0 , HS + dx/10.0  ,  0 , dx/2.0 ,rhoS, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==3)
		{
			dom.Particles[a]->n			= n;
			dom.Particles[a]->k			= k;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->RhoF		= rhoF;
			dom.Particles[a]->IsFree	= false;
		}
	}

	dom.Solve(/*tf*/50000.0,/*dt*/T,/*dtOut*/0.1,"test06",10000);
	return 0;
}
MECHSYS_CATCH
