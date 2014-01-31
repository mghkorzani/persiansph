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

#include <Source/Domain.h>


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
        SPH::Domain dom;
        dom.Gravity = 0.0,-9.81,0.0;
        dom.Dimension = 3;
        dom.Alpha = 1.0;
        dom.Beta = 1.0;
        dom.MaxVel = sqrt(2*9.81*0.00985);
        dom.AutoSaveInt = 1.0;
        dom.Cellfac= 2;
        size_t Nproc = 8;

        dom.AddBoxLength(1,Vec3_t ( 1.001, 0.999 ,0.0  ),0.044,  0,  0, 22 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 0.998, 0.998 ,0.0  ),0.046,  0,  0, 23 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(2,Vec3_t ( 1.045, 0.999 ,0.0  ),0.012,  0,  0,  6 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(2,Vec3_t ( 1.044, 0.998 ,0.0  ),0.014,  0,  0,  7 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 1.057, 0.999 ,0.0  ),0.044,  0,  0, 22 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 1.058, 0.998 ,0.0  ),0.046,  0,  0, 23 ,  1 ,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 0.999, 0.999 ,0.0  ), 0.0 ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 0.998, 1.0   ,0.0  ), 0.0 ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 1.101, 0.999 ,0.0  ), 0.0 ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(1,Vec3_t ( 1.102, 1.0   ,0.0  ), 0.0 ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
        dom.AddRandomBox(3,Vec3_t ( 1.0  , 1.0   ,0.0  ),0.1,0.200,  0, 50 , 100,  1, 0, 1000, 0.002);

//        dom.CellInitiate();
//        dom.ListGenerate();
//        dom.InitiateInteractions();
//        dom.WriteXDMF("test02");
        dom.Solve(/*tf*/0.001,/*dt*/0.0001,/*dtOut*/0.01,"test02",Nproc);
        return 0;
}
MECHSYS_CATCH
