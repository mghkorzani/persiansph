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

void UserAcc(SPH::Domain & domi)
{
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->ID == 3)
			domi.Particles[i]->a = Vec3_t(1.0,0.0,0.0);
		if (domi.Particles[i]->ID == 2)
			domi.Particles[i]->a = Vec3_t(-1.0,0.0,0.0);
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Nproc		= 8;
    	dom.PresEq		= 0;
    	dom.KernelType	= 0;
    	dom.Shepard		= false;

    	dom.TI			= 0.3;
    	dom.Alpha		= 1.0;
    	dom.XSPH		= 0.5;

        double dx,h,rho,K,G;
    	double H,L,n;

    	H	= 0.01;
    	L	= 0.03;
    	n	= 60.0;

    	rho	= 1000.0;
    	K	= 1.0e6;
    	G	= 3.0e5;
    	dx	= H / n;
    	h	= dx*1.5;
        dom.Cs			= sqrt(K/rho);
    	dom.InitialDist	= dx;

        double timestep;
        timestep = (0.1*h/(dom.Cs));
        cout<<timestep<<endl;
        cout<<dom.Cs<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;

     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-L/20.0 , -H/2.0 , 0.0 ), L + L/10.0 + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

     	double x;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G = G;
    		dom.Particles[a]->Material = 2;
    		dom.Particles[a]->Fail = 1;
    		dom.Particles[a]->Sigmay = 3000.0;
    		x = dom.Particles[a]->x(0);
    		if (x<-L/2.0)
    			dom.Particles[a]->ID=2;
    		if (x>L/2.0)
    			dom.Particles[a]->ID=3;
    	}

//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.0005,"test06",5000);
        return 0;
}
MECHSYS_CATCH
