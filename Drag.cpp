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
        dom.Gravity		= 0.4,0.0,0.0;
        dom.Dimension	= 2;
        dom.Alpha		= 0.0;
        dom.Beta		= 0.0;
        dom.MaxVel		= 0.06;
        dom.AutoSaveInt	= 1.0;
        dom.Cellfac		= 2;
        size_t Nproc	= 8;
        dom.Periodic	= true;
        dom.MU			= 0.001;
        dom.PressureBoundary= false;
        dom.XSPH			= 0.5;
        dom.ConstVelPeriodic= 0.0;

        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 1.001801875 , 0.0 ), 0.0034 , 0 , 0 , 906 , 1 , 1 , 0.0000000140625 , 1000 , 0.000004125 , true);
        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 0.999999125 , 0.0 ), 0.0034 , 0 , 0 , 906 , 1 , 1 , 0.0000000140625 , 1000 , 0.000004125 , true);

        double xa,ya,xb,yb,yc;

        dom.AddRandomBox(3 ,Vec3_t ( 1.000 , 1.000 , 0.0 ), 0.003375 , 0.0018 ,  0 , 900 , 480 , 1 , 0.0000000140625 , 1000, 0.000004125);

        for (size_t a=0; a<dom.Particles.Size(); a++)
        {
        	xb=dom.Particles[a]->x(0);
        	yb=dom.Particles[a]->x(1);
        	if ((xb-1.001125)*(xb-1.001125)+(yb-1.0009)*(yb-1.0009)<=0.000000013)
        		{
        		dom.Particles[a]->ID=4;
        		dom.Particles[a]->IsFree=false;
        		}
        }
        dom.DelParticles(4);

        xa=1.0010125;
        ya=sqrt (0.00000001265625-(xa-1.001125)*(xa-1.001125))+1.0009;

   		while (xa<=1.0012375)
        {
   			dom.AddSingleParticle(4,Vec3_t ( xa , ya , 0.0 ), 0.0000000140625 , 1000 , 0.000004125 , true);

   			xb= xa+0.00000000001;
			yb=sqrt (0.00000001265625-(xb-1.001125)*(xb-1.001125))+1.0009;
			std::cout<<yb<<std::endl;
			while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.00000375)
        		{
					xb=xb+0.00000000001;
					yb=sqrt (0.00000001265625-(xb-1.001125)*(xb-1.001125))+1.0009;
       			}
			xa=xb;
			ya=yb;
			yc=(-sqrt (0.00000001265625-(xa-1.001125)*(xa-1.001125))+1.0009);

			if (yc>0.0) dom.AddSingleParticle(4,Vec3_t ( xa , yc , 0.0 ), 0.0000000140625 , 1000 , 0.000004125 , true);
        }
//   		dom.AddSingleParticle(4,Vec3_t ( 1.0012375 , 1.0009 , 0.0 ), 0.0000000140625 , 1000 , 0.000004125 , true);
//
        dom.WriteXDMF("maz");
//
//        dom.Solve(/*tf*/10.0,/*dt*/0.00001,/*dtOut*/0.005,"test07",Nproc);
        return 0;
}
MECHSYS_CATCH
