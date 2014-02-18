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
        dom.Gravity		= 0.01,0.0,0.0;
        dom.Dimension	= 3;
        dom.Alpha		= 0.3;
        dom.Beta		= 0.0;
        dom.MaxVel		= 0.05;
        dom.AutoSaveInt	= 1.0;
        dom.Cellfac		= 2;
        size_t Nproc	= 8;
        dom.Periodic	= true;
        dom.MU			= 0.00089;
        dom.PressureBoundary= false;
        dom.XSPH			= 0.5;

        dom.AddBoxLength(1,Vec3_t ( 1.000 , 1.050 ,0.0  ), 0.2 ,  0,  0,100 ,  1 ,  1, 0, 1000, 0.002, true);
//        dom.AddBoxLength(1,Vec3_t ( 1.001 , 1.051 ,0.0  ),0.198,  0,  0,99  ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 1.000 , 0.999 ,0.0  ), 0.2 ,  0,  0,100 ,  1 ,  1, 0, 1000, 0.002, true);
//        dom.AddBoxLength(1,Vec3_t ( 1.001 , 0.998 ,0.0  ),0.198,  0,  0,99  ,  1 ,  1, 0, 1000, 0.002, true);
//        dom.AddBoxLength(1,Vec3_t ( 1.050, 1.045 ,0.0  ),  0 ,0.004,  0, 1  ,  2 ,  1, 0, 1000, 0.002, true);

//        double maz=1.000;
//        for (size_t i=0; i<100; i+=2)
//        {
//            dom.AddBoxLength(1,Vec3_t ( 1.0 , maz ,0.0  ), 0.2 ,  0,  0,100 ,  1 ,  1, 0.00000201, 1000, 0.0026, false);
//            dom.AddBoxLength(1,Vec3_t ( 1.001 , maz+0.001 ,0.0  ), 0.198 ,  0,  0,99 ,  1 ,  1, 0.00000201, 1000, 0.0026, false);
//             maz+=0.002;
//        }
//        dom.AddRandomBox(3,Vec3_t ( 1.0  , 0.9995   ,0.0  ), 0.2 ,0.05,  0,200 ,100 ,  1, 0.00000025, 1000, 0.002);

        dom.AddRandomBox(3,Vec3_t ( 1.0  , 0.9995   ,0.0  ), 0.2 ,0.05,  0,100 ,25 ,  1, 0, 1000, 0.002);
//        dom.AddRandomBoxFixed(3,Vec3_t ( 1.0  , 0.9995   ,0.0  ), 0.2 ,0.05,  0,200 ,100 ,  1, 0.00000025, 1000, 0.0016, 0.001);

//        dom.WriteXDMF("test05");
        dom.Solve(/*tf*/10.0,/*dt*/0.0001,/*dtOut*/0.005,"test05",Nproc);
        return 0;
}
MECHSYS_CATCH
