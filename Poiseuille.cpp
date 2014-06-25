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
        dom.Gravity		= 0.2,0.0,0.0;
        dom.Dimension	= 2;
        dom.Cs			= 0.1;
        dom.Periodic	= true;
        dom.MU			= 0.001;

        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 1.0010125 , 0.0 ), 0.00375 , 0 , 0 , 150 , 1 , 1 , 0.000000625 , 1000 , 0.0000275 , true);
        dom.AddBoxLength(1 ,Vec3_t ( 1.0000 , 0.9999875 , 0.0 ), 0.00375 , 0 , 0 , 150 , 1 , 1 , 0.000000625 , 1000 , 0.0000275 , true);

        dom.AddBoxLength(1 ,Vec3_t ( 1.0005 , 1.0005 , 0.0 ), 0 , 0.000050 , 0 , 1 , 2 , 1 , 0.000000625 , 1000 , 0.0000275 , true);

        dom.AddRandomBox(3 ,Vec3_t ( 1.000 , 1.000 , 0.0 ), 0.003725 , 0.001 ,  0 , 149 , 40 , 1 , 0.000000625 , 1000, 0.0000275);

        dom.Solve(/*tf*/10.0,/*dt*/0.00001,/*dtOut*/0.002,"test05");
        return 0;
}
MECHSYS_CATCH
