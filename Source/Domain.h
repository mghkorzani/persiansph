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

#ifndef MECHSYS_SPH_DOMAIN_H
#define MECHSYS_SPH_DOMAIN_H

// Std Lib
#include <stdio.h>			/// for NULL
#include <algorithm>		/// for min,max


#include <Source/Interaction.h>

// HDF File Output
#include <hdf5.h>
#include <hdf5_hl.h>

#include <mechsys/util/stopwatch.h>
#include <mechsys/util/string.h>

namespace SPH {



class Domain
{
public:

    // Constructor
    Domain();  ///< Constructor with a vector containing the number of divisions per length, and Xmin Xmax defining the limits of the rectangular domain to be plotted

    // Destructor
    ~Domain ();

    // Methods
    void AddBox                (int tag, Vec3_t const & x, size_t nx, size_t ny, size_t nz,
    		                  double R, double Mass, double Density, double h, bool Fixed);                                    ///< Add a box of particles (should specify radius of particles)
    void AddBoxLength         (int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    		                  double Mass, double Density, double h, bool Fixed);                                              ///< Add a box of particles with length (calculate radius of particles)
    void AddRandomBox         (int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    		                  double Mass, double Density, double h, size_t RandomSeed=100);                                    ///< Add box of random positioned particles (calculate radius of particles)

//    void StartAcceleration   (Vec3_t const & a = Vec3_t(0.0,0.0,0.0));                                                         ///< Add a fixed acceleration
//    void ComputeAcceleration (double dt);                                                                                     ///< Compute the acceleration due to the other particles
//    void Move                  (double dt);                                                                                     ///< Compute the acceleration due to the other particles

    void InitiateInteractions();                                                                                               ///< Reset the interaction array
    void UpdateInteractions ();																								  ///< Update possible interaction list
    void DelParticles        (int const & Tags);														     				      ///< Delete particles by tags
    void Solve                (double tf, double dt, double dtOut, char const * TheFileKey, size_t Nproc);                   ///< The solving function
    void CellInitiate        ();                                                         										  ///< Size of the domain as a rectangular, make cell and HOC
    void CellReset 	       ();                                                         										  ///< Reset HOC to initial value of -1
    void ListGenerate	       ();                                                         										  ///< Generate linked-list
    void ListandInteractionUpdate();                                                         										  ///< Checks if linked-list needs to be updated
    void NeighbourSearch	   (int q1, int q2, int q3);																		  ///< Find neighbour particle and make interaction for cell(q1,q2,q3)

    void StartAcceleration  (Vec3_t const & a);
    void ComputeAcceleration(double dt);
    void Move 				   (double dt);

    void WriteXDMF           (char const * FileKey);                                                                           ///< Save a XDMF file for visualization
    void Save                 (char const * FileKey);                                                         				  ///< Save the domain form a file
    void Load                 (char const * FileKey);                                                         				  ///< Load the domain form a file

    // Data
    Vec3_t                  Gravity;        ///< Gravity acceleration
    Array <Particle*>       Particles;      ///< Array of particles
    Array <Interaction*>    Interactions;   ///< Array of interactions
    Array <Interaction*>    PInteractions;  ///< Array of possible interactions
    double                 Time;           ///< The simulation Time
    size_t                  idx_out;        ///< Index for output purposes
    double 					Dimension;      ///< Dimension of the problem
    double 					Alpha;
    double					Beta;
    double					MaxVel;
    double					AutoSaveInt;	///< Automatic save interval
    Vec3_t                  TRPR;           ///< Top right-hand point at rear of the domain as a Rectangular
    Vec3_t                  BLPF;           ///< Bottom left-hand point At front side of the domain as a Rectangular
    Vec3_t                  CellSize;       ///< Calculated cell size according to (cell size >= 2h)
    int		                CellNo[3];      ///< No. of cells for linked list
    int					*** HOC;			///< Array of (Head of Chain) for each cell
    int					 ** ExInteract;		///< Array to save existing interaction No. of "Interactions" to use in "PInteractions"
    double 					Cellfac;		///< Factor which should be multiplied by h to change the size of cells (min 2)
};
/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Constructor
inline Domain::Domain ()
{
    Time    = 0.0;
    Gravity = 0.0;
    idx_out = 0;
}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactions.Size(); ++i) if (Interactions[i]!=NULL) delete Interactions[i];
}


// Methods
inline void Domain::AddBox(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double R, double Mass, double Density, double h, bool Fixed)
{
    std::cout << "\n--------------Generating particles by AddBox with defined radius-----------------------------------" << std::endl;

    for (size_t i = 0;i<nx;i++)
    for (size_t j = 0;j<ny;j++)
    for (size_t k = 0;k<nz;k++)
    {
        Vec3_t x(2*i*R,2*j*R,2*k*R);
        x +=V;
        if (Mass == 0.0)
        	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,Fixed));
        else
        	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,R,h,Fixed));
     }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, bool Fixed)
{
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length-----------------------------" << std::endl;

	double R = (std::min(Lx/(2*nx),Ly/(2*ny))>0) ? std::min(Lx/(2*nx),Ly/(2*ny)) : std::max(Lx/(2*nx),Ly/(2*ny));
    R = (std::min(R,Lz/(2*nz))>0) ? std::min(R,Lz/(2*nz)) : std::max(R,Lz/(2*nz));

	for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = V(0)+i*Lx/nx;
                double y = V(1)+j*Ly/ny;
                double z = V(2)+k*Lz/nz;
                if (Mass == 0.0)
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,Fixed));
                else
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,R,h,Fixed));
            }
        }
    }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::AddRandomBox(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, size_t RandomSeed)
{
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating random packing of particles by AddRandomBox-------------------------------" << std::endl;

    double R = (std::min(Lx/(2*nx),Ly/(2*ny))>0) ? std::min(Lx/(2*nx),Ly/(2*ny)) : std::max(Lx/(2*nx),Ly/(2*ny));
    R = (std::min(R,Lz/(2*nz))>0) ? std::min(R,Lz/(2*nz)) : std::max(R,Lz/(2*nz));

	std::cout << "\n Radius of Particle = " << R << std::endl;

    double qin = 0.95;
    srand(RandomSeed);
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = V(0)+(i+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lx/nx;
                double y = V(1)+(j+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Ly/ny;
                double z = V(2)+(k+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lz/nz;
                if (Mass == 0.0)
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,false));
                else
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,R,h,false));
            }
        }
    }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::CellInitiate ()
{
	// Calculate Domain Size
	double h=0.0;
	BLPF = Particles[0]->x;
	TRPR = Particles[0]->x;

	for (size_t i=0; i<Particles.Size(); i++)
    {
        if (Particles[i]->x(0) > TRPR(0)) TRPR(0) = Particles[i]->x(0);
        if (Particles[i]->x(1) > TRPR(1)) TRPR(1) = Particles[i]->x(1);
        if (Particles[i]->x(2) > TRPR(2)) TRPR(2) = Particles[i]->x(2);

        if (Particles[i]->x(0) < BLPF(0)) BLPF(0) = Particles[i]->x(0);
        if (Particles[i]->x(1) < BLPF(1)) BLPF(1) = Particles[i]->x(1);
        if (Particles[i]->x(2) < BLPF(2)) BLPF(2) = Particles[i]->x(2);

        if (Particles[i]->h > h) h=Particles[i]->h;
    }
    if ((BLPF(0) < 0.0) | (BLPF(0) < 0.0) | (BLPF(0) < 0.0))
    	{
    	std::cout << "\nProblem to allocate Cells !!!!!!!!!!!!!" << std::endl;
    	std::cout << "Particle with minus coordinate" << std::endl;
    	abort();
    	}

    // Calculate Cells Properties
    CellNo[0] = (int) (floor((TRPR(0)-BLPF(0))/(Cellfac*h)));
    CellNo[1] = (int) (floor((TRPR(1)-BLPF(1))/(Cellfac*h)));
    CellNo[2] = (int) (floor((TRPR(2)-BLPF(2))/(Cellfac*h)));
    CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);

    if (CellNo[2]==0) CellNo[2]=1;
//    std::cout << "h= " << h << std::endl;
//    std::cout << "cell size = " << CellSize << std::endl;
//    std::cout << "BLPF x= " << BLPF(0) << std::endl;
//    std::cout << "BLPF y= " << BLPF(1) << std::endl;
//    std::cout << "BLPF z= " << BLPF(2) << std::endl;
//    std::cout << "TRPR x= " << TRPR(0) << std::endl;
//    std::cout << "TRPR y= " << TRPR(1) << std::endl;
//    std::cout << "TRPR z= " << TRPR(2) << std::endl;
//
    std::cout << "No of Cells in X Direction = " << CellNo[0] << std::endl;
    std::cout << "No of Cells in Y Direction = " << CellNo[1] << std::endl;
    std::cout << "No of Cells in Z Direction = " << CellNo[2] << std::endl;

    // Initiate Head of Chain array for Linked-List
    HOC = new int**[(int) CellNo[0]];
    for(int i =0; i<CellNo[0]; i++){
       HOC[i] = new int*[CellNo[1]];
       for(int j =0; j<CellNo[1]; j++){
           HOC[i][j] = new int[CellNo[2]];
           for(int k = 0; k<CellNo[2];k++){
              HOC[i][j][k] = -1;
//              std::cout << i <<j << k << std::endl;

           }
       }
    }

}

inline void Domain::CellReset ()
{
    for(int i =0; i<CellNo[0]; i++){
       for(int j =0; j<CellNo[1]; j++){
           for(int k = 0; k<CellNo[2];k++){
              HOC[i][j][k] = -1;
           }
       }
    }

}

inline void Domain::ListGenerate ()
{
	int i,j,k;
	k=0;
	double temp;
    for (size_t a=0; a<Particles.Size(); a++)
    {
        i= (int) (Particles[a]->x(0) - BLPF(0)) / CellSize(0);
        j= (int) (Particles[a]->x(1) - BLPF(1)) / CellSize(1);
//        k= (int) Particles[a]->x(2) / CellSize(2);
        k=0;

        temp = HOC[i][j][k];
        HOC[i][j][k] = a;
        Particles[a]->LL = temp;
        Particles[a]->CC[0] = i;
        Particles[a]->CC[1] = j;
        Particles[a]->CC[2] = k;
//        std::cout << HOC[i][j][k] <<"  "<< Particles[a]->LL <<std::endl;
    }

}

inline void Domain::ListandInteractionUpdate()
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
	if (Particles[i]->CellUpdate(CellSize,BLPF))
		{
		CellReset();
	    for (size_t a=0; a<Particles.Size(); a++)
	    	{
	        Particles[a]->LL = -1;
	    	}
	    ListGenerate();
	    UpdateInteractions();
	    break;
		}
    }
}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->a = a;
        Particles[i]->dDensity = 0.0;
    }
}

inline void Domain::ComputeAcceleration (double dt)
{
	#pragma omp parallel for
	for (size_t i=0; i<PInteractions.Size(); i++) PInteractions[i]->CalcForce(dt);
}

inline void Domain::Move (double dt)
{
   	#pragma omp parallel for
	for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt);
}

inline void Domain::DelParticles (int const & Tags)
{
    Array<int> idxs; // indices to be deleted
    for (size_t i=0; i<Particles.Size(); ++i)
    	{
        if (Particles[i]->ID==Tags) idxs.Push(i);
    	}
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particle to be deleted");
    Particles.DelItems (idxs);
    std::cout << "\n" << "Particles with Tag No. " << Tags << " has been deleted" << std::endl;

}

inline void Domain::NeighbourSearch(int q1, int q2, int q3)
{
	int temp1, temp2;
	temp1 = HOC[q1][q2][q3];

	while (temp1 != -1)
	{

		// The current cell  => self interactions
		temp2 = Particles[temp1]->LL;
		while (temp2 != -1)
		{
			if (ExInteract [temp1][temp2] == -1)
			{
				Interactions.Push(new Interaction(Particles[temp1],Particles[temp2],Dimension,Alpha,Beta,MaxVel));
				ExInteract [temp1][temp2] = Interactions.size() - 1;
				ExInteract [temp2][temp1] = Interactions.size() - 1;
				PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
			}
			else
			{
				PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
				}
			temp2 = Particles[temp2]->LL;
		}


		// (q1 + 1, q2 , q3)
		if (q1+1< CellNo[0])
		{
			temp2 = HOC[q1+1][q2][q3];
			while (temp2 != -1)
			{
				if (ExInteract [temp1][temp2] == -1)
				{
					Interactions.Push(new Interaction(Particles[temp1],Particles[temp2],Dimension,Alpha,Beta,MaxVel));
					ExInteract [temp1][temp2] = Interactions.size() - 1;
					ExInteract [temp2][temp1] = Interactions.size() - 1;
					PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
				}
				else
				{
					PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
				}
				temp2 = Particles[temp2]->LL;
			}
		}

		// (q1 + a, q2 + 1, q3) & a[-1,1]
		if (q2+1< CellNo[1])
		{
			for (int i = q1-1; i <= q1+1; i++)
			{
				if (i<CellNo[0] && i>=0)
				{
					temp2 = HOC[i][q2+1][q3];
					while (temp2 != -1)
					{
						if (ExInteract [temp1][temp2] == -1)
						{
							Interactions.Push(new Interaction(Particles[temp1],Particles[temp2],Dimension,Alpha,Beta,MaxVel));
							ExInteract [temp1][temp2] = Interactions.size() - 1;
							ExInteract [temp2][temp1] = Interactions.size() - 1;
							PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
						}
						else
						{
							PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
						}
					temp2 = Particles[temp2]->LL;
					}
				}
			}
		}

		// (q1 + a, q2 + b, q3 + 1) & a,b[-1,1] => all 9 cells above the current cell
		if (q3+1< CellNo[2])
		{
			for (int j=q2-1; j<=q2+1; j++)
			{
				for (int i=q1-1; i<=q1+1; i++)
				{
					if (i<CellNo[0] && i>=0 && j<CellNo[1] && j>=0)
					{
						temp2 = HOC[i][j][q3+1];
						while (temp2 != -1)
						{
							if (ExInteract [temp1][temp2] == -1)
							{
								Interactions.Push(new Interaction(Particles[temp1],Particles[temp2],Dimension,Alpha,Beta,MaxVel));
								ExInteract [temp1][temp2] = Interactions.size() - 1;
								ExInteract [temp2][temp1] = Interactions.size() - 1;
								PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
							}
							else
							{
								PInteractions.Push(Interactions[ ExInteract [temp1][temp2] ]);
							}
							temp2 = Particles[temp2]->LL;
						}
					}
				}
			}
		}

		temp1 = Particles[temp1]->LL;
	}
}

inline void Domain::InitiateInteractions()
{
	// Initiate ExInteract
	ExInteract = new int*[Particles.Size()];
    for(size_t i =0; i<Particles.Size(); i++)
    	{
    	ExInteract[i] = new int[Particles.Size()];
    	for(size_t j =0; j<Particles.Size(); j++)
    		{
    		ExInteract[i][j] = -1;
    		}
    	}

	// Delete old interactions
    for (size_t i=0; i<Interactions.Size(); ++i)
    {
        if (Interactions[i]!=NULL) delete Interactions[i];
    }
    PInteractions.Resize(0);

    for (int q3=0; q3<CellNo[2]; q3++)
    {
        for (int q2=0; q2<CellNo[1]; q2++)
        {
        	for (int q1=0; q1<CellNo[0]; q1++)
            {
            	if (HOC[q1][q2][q3]==-1) continue;
            	else  NeighbourSearch(q1,q2,q3);
            }
        }
    }
}

inline void Domain::UpdateInteractions()
{
    PInteractions.Resize(0);

    for (int q3=0; q3<CellNo[2]; q3++)
    {
        for (int q2=0; q2<CellNo[1]; q2++)
        {
            for (int q1=0; q1<CellNo[0]; q1++)
            {
            	if (HOC[q1][q2][q3]==-1) continue;
            	else  NeighbourSearch(q1,q2,q3);
            }

        }

    }
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t Nproc)
{
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Solving------------------------------------------------------------------------------" << std::endl;

    idx_out = 1;
    double tout = Time;

    size_t save_out = 1;
    double sout = AutoSaveInt;

    CellInitiate();
    ListGenerate();
    InitiateInteractions();

    while (Time<tf)
    {

    	StartAcceleration(Gravity);
    	ComputeAcceleration(dt);
    	Move(dt);


        // output
        if (Time>=tout)
        {
            if (TheFileKey!=NULL)
            	{
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                WriteXDMF    (fn.CStr());
                std::cout << "\n" << "Output No. " << idx_out << " at " << Time << " has been generated" << std::endl;
            	}
            idx_out++;
            tout += dtOut;
        }

        // Auto Save
       if (AutoSaveInt>0)
       {
    	   if (Time>=sout)
    	   {
    		   if (TheFileKey!=NULL)
    		   {
               String fn;
               fn.Printf    ("Auto Save_%s_%04d", TheFileKey, save_out);
               Save   		(fn.CStr());
               std::cout << "\n" << "Auto Save No. " << save_out << " at " << Time << " has been generated" << std::endl;
    		   }
    		save_out++;
    		sout += AutoSaveInt;
    	   }
       }

       Time += dt;

       ListandInteractionUpdate();

    }
}

inline void Domain::WriteXDMF (char const * FileKey)
{

	String fn(FileKey);
    fn.append(".h5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    float * Posvec   = new float[3*Particles.Size()];
    float * Velvec   = new float[3*Particles.Size()];
    float * Pressure = new float[  Particles.Size()];
    float * Radius   = new float[  Particles.Size()];
    int   * Tag      = new int  [  Particles.Size()];

    for (size_t i=0;i<Particles.Size();i++)
    {
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        Pressure[i    ] = float(Particles[i]->Pressure);
        Radius  [i    ] = float(Particles[i]->R);
        Tag     [i    ] = int  (Particles[i]->ID);
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = Particles.Size();
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);


    delete [] Posvec;
    delete [] Velvec;
    delete [] Pressure;
    delete [] Radius;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);


    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"SPHCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::Save (char const * FileKey)
{
    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");

    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing the number of particles in the domain
    int data[1];
    data[0]=Particles.Size();
    hsize_t dims[1];
    dims[0]=1;
    H5LTmake_dataset_int(file_id,"/NP",1,dims,data);

    for (size_t i=0; i<Particles.Size(); i++)
    {
        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gcreate2(file_id, par.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Storing some scalar variables
        double dat[1];
        dat[0] = Particles[i]->Mass;
        H5LTmake_dataset_double(group_id,"Mass",1,dims,dat);
        dat[0] = Particles[i]->Density;
        H5LTmake_dataset_double(group_id,"Rho",1,dims,dat);
        dat[0] = Particles[i]->hr;
        H5LTmake_dataset_double(group_id,"h",1,dims,dat);
        dat[0] = Particles[i]->R;
        H5LTmake_dataset_double(group_id,"Radius",1,dims,dat);

        int tag[1];
        tag[0] = Particles[i]->ID;
        H5LTmake_dataset_int(group_id,"Tag",1,dims,tag);

        // Change  to integer for fixity
        int IsFree[1];
        if (Particles[i]->IsFree == true) IsFree[0]=1;
        		else IsFree[0]=0;
        H5LTmake_dataset_int(group_id,"IsFree",1,dims,IsFree);

        // Storing vectorial variables
        double cd[3];
        hsize_t dd[1];
        dd[0] = 3;

        cd[0]=Particles[i]->x(0);
        cd[1]=Particles[i]->x(1);
        cd[2]=Particles[i]->x(2);
        H5LTmake_dataset_double(group_id,"x",1,dd,cd);

        cd[0]=Particles[i]->v(0);
        cd[1]=Particles[i]->v(1);
        cd[2]=Particles[i]->v(2);
        H5LTmake_dataset_double(group_id,"v",1,dd,cd);
    }

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
}

inline void Domain::Load (char const * FileKey)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"/NP",data);
    size_t NP = data[0];

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gopen2(file_id, par.CStr(),H5P_DEFAULT);

        Particles.Push(new Particle(0,Vec3_t(0,0,0),Vec3_t(0,0,0),0,0,0,0,false));

        // Loading vectorial variables
        double cd[3];
        H5LTread_dataset_double(group_id,"x",cd);
        Particles[Particles.Size()-1]->x = Vec3_t(cd[0],cd[1],cd[2]);
        Particles[Particles.Size()-1]->xb = Vec3_t(cd[0],cd[1],cd[2]);			// Because of the constructor in Particle

//        H5LTread_dataset_double(group_id,"v",cd);
//        Particles[Particles.Size()-1]->v = Vec3_t(cd[0],cd[1],cd[2]);

        // Loading the scalar quantities of the particle
        double dat[1];
        H5LTread_dataset_double(group_id,"Mass",dat);
        Particles[Particles.Size()-1]->Mass = dat[0];

        H5LTread_dataset_double(group_id,"Rho",dat);
        Particles[Particles.Size()-1]->Density = dat[0];
        Particles[Particles.Size()-1]->RefDensity = dat[0];			// Because of the constructor in Particle
        Particles[Particles.Size()-1]->Densityb = dat[0];			// Because of the constructor in Particle

        H5LTread_dataset_double(group_id,"h",dat);
        Particles[Particles.Size()-1]->hr = dat[0];
        Particles[Particles.Size()-1]->h = dat[0];			// Because of the constructor in Particle

        H5LTread_dataset_double(group_id,"Radius",dat);
        Particles[Particles.Size()-1]->R = dat[0];

        int datint[1];
        H5LTread_dataset_int(group_id,"Tag",datint);
        Particles[Particles.Size()-1]->ID = datint[0];

        H5LTread_dataset_int(group_id,"IsFree",datint);
        if (datint[0] == 1) Particles[Particles.Size()-1]->IsFree=true;
        		else Particles[Particles.Size()-1]->IsFree=false;
    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

}; // namespace SPH

#endif // MECHSYS_SPH_DOMAIN_H
