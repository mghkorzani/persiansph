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
        dom.Gravity		= 0.0002,0.0,0.0;
        dom.Dimension	= 2;
        dom.Alpha		= 0.0;
        dom.Beta		= 0.0;
        dom.MaxVel		= 0.000025;
        dom.AutoSaveInt	= 1.0;
        dom.Cellfac		= 2;
        size_t Nproc	= 8;
        dom.Periodic	= true;
        dom.MU			= 0.001;
        dom.PressureBoundary= false;
        dom.XSPH			= 0.0;
        dom.ConstVelPeriodic= 0.0;

        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 1.0010125 , 0.0 ), 0.000525 , 0 , 0 , 21 , 1 , 1 , 0.000000625 , 1000 , 0.0000275 , true);
        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 0.9999875 , 0.0 ), 0.000525 , 0 , 0 , 21 , 1 , 1 , 0.000000625 , 1000 , 0.0000275 , true);

        dom.AddRandomBox(3 ,Vec3_t ( 1.000 , 1.000 , 0.0 ), 0.0005 , 0.001 ,  0 , 20 , 40 , 1 , 0.000000625 , 1000, 0.0000275);

        dom.Solve(/*tf*/10.0,/*dt*/0.00001,/*dtOut*/0.001,"test07",Nproc);
        return 0;
}
MECHSYS_CATCH
