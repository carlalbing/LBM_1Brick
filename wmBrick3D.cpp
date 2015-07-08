#include <mpi.h>

#include <iostream>
//#include "OpenChannel3D.h"
#include "WallMountedBrick.h"
#include "workArounds.h"

using namespace std;

int main(int argc, char * argv[])
{
  const char * paramFN;
  const char * obstFN;
  int rank, size;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  paramFN = "params.lbm";
  obstFN = "obst_file.lbm";
  if (argc >= 2) { paramFN = argv[1]; } 
  if (argc >= 3) { obstFN = argv[2]; }

  WallMountedBrick pp(rank, size, paramFN, obstFN);

  if(rank==0){
	  cout<< "Problem initialized!" << endl;
  }

  int* inl = pp.inl;
  int* onl = pp.onl;
  int* snl = pp.snl;
  float* u_bc = pp.u_bc;
  int nnodes = pp.nnodes;

  dummyUse(inl, onl, snl, u_bc, nnodes);

  #pragma acc  data \
      copyin(inl[0:nnodes], onl[0:nnodes], snl[0:nnodes], u_bc[0:nnodes])
  {
    // write initial data
    // (data processing script will be expecting it)
    pp.write_data(MPI_COMM_WORLD,true);

    double time_start, time_end, ex_time, LPU_sec, gNumLP;
    time_start = MPI_Wtime();
    for(int ts = 0; ts<pp.Num_ts;ts++){
	    // say something comforting about the problem progress
	    if((ts+1)%(pp.ts_rep_freq)==0){
		    if(rank==0){
			    cout << "Executing time step number " << ts+1 << endl;
		    }

	    }
	    pp.take_lbm_timestep(ts%2==0,MPI_COMM_WORLD); // weird function call sig.

	     if((ts+1)%(pp.plot_freq)==0){
	     	  // write data at requested intervals.
	     	  pp.write_data(MPI_COMM_WORLD,ts%2==0);
	     }
    }
    time_end = MPI_Wtime();

    if(rank==0){
	    ex_time = time_end - time_start;
	    gNumLP = pp.Nx*pp.Ny*pp.Nz;
	    LPU_sec = ((double)gNumLP*(double)pp.Num_ts)/ex_time;
	    cout << "Estimiated LPU/sec = " << LPU_sec << endl;
    }
  }
  MPI_Finalize();
  return 0;
}
