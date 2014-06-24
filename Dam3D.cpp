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

	dom.Dimension	= 3;
	dom.Gravity     = 0.0,-9.81,0.0;
	dom.Cs			= 1.5;
	dom.Alpha		= 0.05;
	dom.P0			= 0.0;
	dom.PresEq		= 0;
	dom.Nproc		= 1;
	dom.TI			= 0.1;

	dom.AddBoxLength(1 ,Vec3_t ( 1.0 , 1.0 , 1.0                ), 0.0045 , 0.00225 , 0               , 100 , 50 , 1 , 9.1125e-11 , 1000 , 5.00e-5 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.0 , 1.0-2.25e-5 , 1.0-2.25e-5        ), 0.0045 , 0.00225 , 0               , 100 , 50 , 1 , 9.1125e-11 , 1000 , 5.00e-5 , true);

	dom.AddBoxLength(2 ,Vec3_t ( 1.0 , 1.0 , 1.001125           ), 0.0045 , 0.00225 , 0               , 100 , 50 , 1 , 9.1125e-11 , 1000 , 5.00e-5 , true);
	dom.AddBoxLength(2 ,Vec3_t ( 1.0 , 1.0-2.25e-5 , 1.001125+2.25e-5   ), 0.0045 , 0.00225 , 0               , 100 , 50 , 1 , 9.1125e-11 , 1000 , 5.00e-5 , true);

	dom.AddBoxLength(3 ,Vec3_t ( 1.0 , 1.0-2.25e-5 , 1.0        ), 0.0045 , 0       , 0.001125+4.50e-5, 100 , 1 , 26 , 9.1125e-11 , 1000 , 5.00e-5 , true);
	dom.AddBoxLength(3 ,Vec3_t ( 1.0-2.25e-5 , 1.0-4.50e-5 , 1.0      ), 0.004545 , 0       , 0.001125+4.50e-5, 101 , 1 , 26 , 9.1125e-11 , 1000 , 5.00e-5 , true);

	dom.AddBoxLength(4 ,Vec3_t ( 1.0-2.25e-5 , 1.0 , 1.0), 0      , 0.00225 , 0.001125+4.50e-5 , 1 , 50 , 26 , 9.1125e-11 , 1000 , 5.00e-5 , true);
	dom.AddBoxLength(4 ,Vec3_t ( 1.0-4.5e-5, 1.0 , 1.0-2.25e-5), 0      , 0.00225 , 0.001125+2*4.50e-5 , 1 , 50 , 27 , 9.1125e-11 , 1000 , 5.00e-5 , true);

	dom.AddBoxLength(5 ,Vec3_t ( 1.0045-2.25e-5 , 1.0-4.50e-5 , 1.0     ), 0      , 0.00225 , 0.001125+4.50e-5 , 1 , 50 , 26 , 9.1125e-11 , 1000 , 5.00e-5 , true);
	dom.AddBoxLength(5 ,Vec3_t ( 1.0045, 1.0-4.50e-5, 1.0-2.25e-5), 0      , 0.00225 , 0.001125+2*4.50e-5 , 1 , 50 , 27 , 9.1125e-11 , 1000 , 5.00e-5 , true);

	dom.AddRandomBox(6 ,Vec3_t ( 1.0000225 , 1.0000225 , 1.0    ), 0.00225, 0.0018,  0.00045         ,50 , 40  , 10, 9.1125e-11 , 1000, 5.00e-5);

//	dom.WriteXDMF("maz");

	dom.Solve(/*tf*/0.005,/*dt*/0.00001,/*dtOut*/0.0005,"test07");
	return 0;
}
MECHSYS_CATCH
