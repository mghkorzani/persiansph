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
	dom.Gravity = 0.0,-9.81,0.0;

	dom.AddBoxLength(-1,Vec3_t ( 0.0 , 0.0 ,0.0  ), 50,  0,  0,  50,  1 ,  1, 0, 1000, 0.5, true);
	dom.AddBoxLength(-1,Vec3_t ( 0.0 , 1.0 ,0.0  ), 0 , 50,  0, 1  ,  50,  1, 0, 1000, 0.5, true);
	dom.AddBoxLength(-1,Vec3_t ( 0.0 , 51.0,0.0  ), 50,  0,  0,  50,  1 ,  1, 0, 1000, 0.5, true);
	dom.AddBoxLength(-1,Vec3_t ( 49.0, 1.0 ,0.0  ), 0 , 50,  0, 1  ,  50,  1, 0, 1000, 0.5, true);
	dom.AddRandomBox(-2,Vec3_t ( 0.5 , 0.5 ,0.0  ), 20, 40,  0, 20 ,  40,  1, 0, 1000, 0.5);
//	dom.WriteXDMF("test02");
	dom.Solve(/*tf*/60.0,/*dt*/0.001,/*dtOut*/0.05,"test02");
	return 0;
}
MECHSYS_CATCH
