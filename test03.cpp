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
        dom.MaxVel = sqrt(2*9.81*0.01);
        dom.AutoSaveInt = 1.0;
        size_t Nproc = 8;

        dom.Load("Auto Save_test02_0002");
        dom.AddBoxLength(4,Vec3_t (-0.049,-0.151,0.0 ), 0.2 ,  0 ,  0, 100,  1,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(4,Vec3_t (-0.052,-0.152,0.0 ),0.206,  0 ,  0, 103,  1,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(5,Vec3_t (-0.051,-0.151,0.0 ), 0.0 ,0.15,  0, 1  , 75,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(5,Vec3_t (-0.052,-0.150,0.0 ), 0.0 ,0.15,  0, 1  , 75,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(6,Vec3_t ( 0.151,-0.151,0.0 ), 0.0 ,0.15,  0, 1  , 75,  1, 0, 1000, 0.002, true);
        dom.AddBoxLength(6,Vec3_t ( 0.152,-0.150,0.0 ), 0.0 ,0.15,  0, 1  , 75,  1, 0, 1000, 0.002, true);
        dom.DelParticles(2);

//        dom.WriteXDMF("test02");

        dom.Solve(/*tf*/5.0,/*dt*/0.0001,/*dtOut*/0.005,"test02",Nproc);
        return 0;
}
MECHSYS_CATCH
