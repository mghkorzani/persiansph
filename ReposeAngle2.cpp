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
        dom.Nproc		= 24;
    	dom.KernelType	= 0;
    	dom.Alpha		= 0.1;
    	dom.Beta		= 0.1;
    	dom.Gravity		= 0.0, -9.81, 0.0;
//    	dom.TI			= 0.5;
//    	dom.TIn			= 2.55;
// 		dom.Scheme		= 1; //Leapfrog
		dom.Scheme		= 0; //Modified Verlet

        double dx,h,rho,K,G,timestep,E,Nu,H,n,L,Phi,Psi,c,Cs;

    	H	= 2.0;
    	n	= 50.0;
    	L	= 4.0;
    	rho	= 1850.0;
    	K	= 1.5e6;
    	Phi	= 25.0;
    	c	= 0.0;
    	Psi	= 0.0;
    	Nu	= 0.3;
    	E	= (3.0*(1.0-2.0*Nu))*K;
    	G	= E/(2.0*(1.0+Nu));
    	dx	= H / n;
    	h	= dx*1.2;
        Cs	= sqrt(E/rho);
        cout<<"Cs = "<<Cs<<endl;
        Cs	= 600.0/6.0; //Very important
    	dom.InitialDist	= dx;
        timestep		= (0.2*h/(Cs));

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"E = "<<E<<endl;
        cout<<"K = "<<K<<endl;
        cout<<"h = "<<h<<endl;

     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx ,  0.0    , 0.0 ), L + 3.0*dx + dx/10.0 , 1.1*H + dx/10.0  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , 0.0 ), 3.0*L + 3.0*dx      , 3.0*dx + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

		double K0 = 1-sin(Phi/180.0*M_PI);

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G			= G;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->K			= K;
    		dom.Particles[a]->Material	= 3;
    		dom.Particles[a]->Fail		= 3;
    		dom.Particles[a]->c			= c;
    		dom.Particles[a]->phi		= Phi/180.0*M_PI;
    		dom.Particles[a]->psi		= Psi/180.0*M_PI;
    		if (dom.Particles[a]->x(1)<0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}
    		if (dom.Particles[a]->x(0)<0.0 && dom.Particles[a]->x(1)>0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}
    		if (dom.Particles[a]->x(1)>H && dom.Particles[a]->ID==1)
    		{
    			dom.Particles[a]->ID		= 6;
    		}
         	for (size_t i=0; i<10; i++)
    			if (dom.Particles[a]->ID == 1)
    			{
    				if (dom.Particles[a]->x(0)>(9.0*dx + i*10.0*dx + dx/2.1) && dom.Particles[a]->x(0)<=(10.0*dx + i*10.0*dx + dx/2.1))
    					dom.Particles[a]->ID = 2;
    			}

         	for (size_t i=0; i<5; i++)
    			if (dom.Particles[a]->ID == 1)
    			{
    				if (dom.Particles[a]->x(1)>(9.0*dx + i*10.0*dx + dx/2.1) && dom.Particles[a]->x(1)<=(10.0*dx + i*10.0*dx + dx/2.1))
    					dom.Particles[a]->ID = 2;
    			}

			if (dom.Particles[a]->ID == 1 || dom.Particles[a]->ID == 2)
			{
				dom.Particles[a]->Sigma(1,1)	= rho * -9.81 *(H-dom.Particles[a]->x(1));
				dom.Particles[a]->Sigma(0,0)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
				dom.Particles[a]->Sigma(2,2)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
			}
    	}
     	dom.DelParticles(6);

    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.02,"test06",250);
        return 0;
}
MECHSYS_CATCH
