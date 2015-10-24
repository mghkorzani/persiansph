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

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Nproc		= 8;
    	dom.PresEq		= 0;
    	dom.KernelType	= 0;
//    	dom.Shepard		= true;

//    	dom.TI			= 0.3;
//    	dom.TIn			= 4.0;
    	dom.Alpha		= 1.0;
//    	dom.Beta		= 1.0;
//    	dom.XSPH		= 0.5;
    	dom.Gravity		= 0.0, -9.81, 0.0;

        double dx,h,rho,K,G;
    	double E,Nu,H,n;


    	H	= 5.0;
    	n	= 50.0;

    	rho	= 2700.0;
    	E	= 150.0e6;
    	Nu	= 0.3;
    	K	= E/(3.0*(1.0-2.0*Nu));
    	G	= E/(2.0*(1.0+Nu));
    	dx	= H / n;
    	h	= dx*1.5;
        dom.Cs			= sqrt(K/rho);
    	dom.InitialDist	= dx;
    	dom.DomMax(1)	= 1.1*H;

        double timestep;
        timestep = (0.05*h/(dom.Cs));

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<dom.Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"K = "<<K<<endl;

//     	dom.AddBoxLength(1 ,Vec3_t ( -H - 3.0*dx     ,  0.0    , 0.0 ), 2.0*H + 6.0*dx + dx/10.0 , H + dx/10.0      ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -H/2.0     ,  0.0    , 0.0 ), 1.0*H + dx/10.0 , H + dx/10.0      ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -1.0*H , -3.0*dx , 0.0 ), 2.0*H + dx/10.0 , 3.0*dx + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G			= G;
    		dom.Particles[a]->K			= K;
    		dom.Particles[a]->Material	= 3;
    		dom.Particles[a]->Fail		= 2;
    		dom.Particles[a]->c			= 0.0;
    		dom.Particles[a]->phi		= 30.0/180.0*M_PI;
//    		if (dom.Particles[a]->x(0)<-H || dom.Particles[a]->x(0)>H)
//    		{
////    			dom.Particles[a]->NoSlip	= true;
//    			dom.Particles[a]->IsFree	= false;
//    			dom.Particles[a]->ID		= 2;
//    		}
    		if (dom.Particles[a]->x(1)<0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 2;
    		}
//    		dom.Particles[a]->ShearStress(0,0) = rho * -9.81 *(H-dom.Particles[a]->x(1));
//    		dom.Particles[a]->ShearStressb(0,0) = rho * -9.81 *(H-dom.Particles[a]->x(1));
//    		dom.Particles[a]->ShearStress(1,1) = rho * -9.81 *(H-dom.Particles[a]->x(1));
//    		dom.Particles[a]->ShearStressb(1,1) = rho * -9.81 *(H-dom.Particles[a]->x(1));
    	}


    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.01,"test06",999);
        return 0;
}
MECHSYS_CATCH
