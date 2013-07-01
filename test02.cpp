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
	dom.MaxVel = sqrt(2*9.81*0.05);

	dom.AddBoxLength(-1,Vec3_t ( 0.001,-0.001 ,0.0  ), 0.1  ,  0,  0, 50 ,  1 ,  1, 0, 1000, 0.002, true);
//	dom.AddBoxLength(-1,Vec3_t (-0.005,-0.003 ,0.0  ), 0.112,  0,  0, 56 ,  1 ,  1, 0, 1000, 0.005, true);
//	dom.AddBoxLength(-1,Vec3_t (-0.005,-0.005 ,0.0  ), 0.112,  0,  0, 56 ,  1 ,  1, 0, 1000, 0.005, true);
	dom.AddBoxLength(-1,Vec3_t (-0.001,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
//	dom.AddBoxLength(-1,Vec3_t (-0.003,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.005, true);
//	dom.AddBoxLength(-1,Vec3_t (-0.005,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.005, true);
	dom.AddBoxLength(-1,Vec3_t ( 0.101,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.002, true);
//	dom.AddBoxLength(-1,Vec3_t ( 0.103,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.005, true);
//	dom.AddBoxLength(-1,Vec3_t ( 0.105,-0.001 ,0.0  ), 0  ,0.3,  0, 1  , 150,  1, 0, 1000, 0.005, true);
	dom.AddRandomBox(-2,Vec3_t ( 0.0  , 0.0   ,0.0  ), 0.04,0.04,  0, 20 , 20,  1, 0, 1000, 0.002);
//	dom.WriteXDMF("test02");
	dom.Solve(/*tf*/3.0,/*dt*/0.0001,/*dtOut*/0.005,"test02");
	return 0;
}
MECHSYS_CATCH
