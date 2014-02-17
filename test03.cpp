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
    SPH::Domain dom;
    dom.Gravity 		= 0.0,-9.81,0.0;
    dom.Dimension 		= 3;
    dom.Alpha 			= 1.0;
    dom.Beta 			= 1.0;
    dom.MaxVel 			= sqrt(2*9.81*0.2);
    dom.AutoSaveInt 	= 1.0;
    dom.Cellfac			= 2;
    size_t Nproc 		= 8;
    dom.Periodic		= false;
    dom.MU				= 0.0;
    dom.PressureBoundary= false;
    dom.XSPH			= 0.0;

        dom.Load("Auto Save_test02_0029");
        dom.AddBoxLength(1,Vec3_t (-0.099,-0.076,0.0 ), 0.3 ,  0 ,  0, 150,  1,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t (-0.102,-0.077,0.0 ),0.306,  0 ,  0, 153,  1,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t (-0.101,-0.076,0.0 ), 0.0 ,0.10,  0, 1  , 50,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t (-0.102,-0.075,0.0 ), 0.0 ,0.10,  0, 1  , 50,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 0.201,-0.076,0.0 ), 0.0 ,0.10,  0, 1  , 50,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 0.202,-0.075,0.0 ), 0.0 ,0.10,  0, 1  , 50,  1, 0, 1000, 0.002, true);
        dom.DelParticles(2);

//        dom.WriteXDMF("test02");

        dom.Solve(/*tf*/15.0,/*dt*/0.00005,/*dtOut*/0.0005,"test03",Nproc);
        return 0;
}
MECHSYS_CATCH
