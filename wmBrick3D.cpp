#include <mpi.h>

#include <iostream>
//#include "OpenChannel3D.h"
#include "WallMountedBrick.h"
#include "workArounds.h"

#ifdef CRAYPAT
#include "pat_api.h"
#endif

using namespace std;

int main(int argc, char * argv[])
{
     #ifdef CRAYPAT
     PAT_record(PAT_STATE_OFF);
     #endif
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
    
    #ifdef CRAYPAT
    PAT_record(PAT_STATE_ON);
    #endif
    
    int* inl = pp.inl;
    int* onl = pp.onl;
    int* snl = pp.snl;
    float* u_bc = pp.u_bc;
    int nnodes = pp.nnodes;
    
    int * Mspeeds = pp.Mspeeds;
    int * Pspeeds = pp.Pspeeds;
    
    int numPspeeds = pp.numPspeeds;
    int numMspeeds = pp.numMspeeds;
    
    int numSpd = pp.numSpd;
    float * fEven = pp.fEven;
    float * fOdd = pp.fOdd;
    
    float* ux_l = pp.ux_l;
    float* uy_l = pp.uy_l;
    float* uz_l = pp.uz_l;
    float* rho_l = pp.rho_l;
    int numEntries = pp.numEntries;
    
    
    dummyUse(inl, onl, snl, u_bc, nnodes, Mspeeds, Pspeeds, numPspeeds, numMspeeds,numSpd,fEven,fOdd, ux_l, uy_l, uz_l, rho_l, numEntries);
    
    #pragma acc data \
        copyin(inl[0:nnodes], onl[0:nnodes], snl[0:nnodes], u_bc[0:nnodes]) \
        copyin(Mspeeds[0:numMspeeds],Pspeeds[0:numPspeeds]) \
        copyin(fEven[0:nnodes*numSpd]) \
        create(fOdd[0:nnodes*numSpd]) \
        create(ux_l[0:numEntries],uy_l[0:numEntries],uz_l[0:numEntries],rho_l[0:numEntries])        
    {
        // write initial data
        // (data processing script will be expecting it)
        pp.write_data_GPU2Buf(true);
        //pp.write_data_Buf2File(MPI_COMM_WORLD,true);
        
        double time_start, time_end, ex_time, LPU_sec, gNumLP;
        time_start = MPI_Wtime();
        int ts;
        for(ts = 0; ts<pp.Num_ts;ts+=pp.plot_freq){
            #pragma omp parallel num_threads(2)
            #pragma omp single
            {
                #pragma omp task
                for(int ts2 = ts; ts2<ts+pp.plot_freq; ts2++){
                    // say something comforting about the problem progress
                    if((ts2+1)%(pp.ts_rep_freq)==0){
                        if(rank==0){
                            cout << "Executing time step number " << ts2+1 << endl;
                        }
                        
                    }
                    
                    #ifdef _OPENACC
                      pp.take_lbm_timestep_acc(ts2%2==0,MPI_COMM_WORLD); // weird function call sig.
                    #else
                      pp.take_lbm_timestep(ts2%2==0,MPI_COMM_WORLD);
                    #endif
                }
                
                #pragma omp task
                if(ts == 0)
                    pp.write_data_Buf2File(MPI_COMM_WORLD, true);
                else
                    pp.write_data_Buf2File(MPI_COMM_WORLD,(ts-1)%2==0);
            }
            pp.write_data_GPU2Buf((ts+pp.plot_freq-1)%2==0);
        }
        // write data at requested intervals.
        pp.write_data_Buf2File(MPI_COMM_WORLD,(ts-1)%2==0);
        time_end = MPI_Wtime();
        
        if(rank==0){
            ex_time = time_end - time_start;
            gNumLP = pp.Nx*pp.Ny*pp.Nz;
            LPU_sec = ((double)gNumLP*(double)pp.Num_ts)/ex_time;
            cout << "Estimiated LPU/sec = " << LPU_sec << endl;
        }
        
    } // end of acc data region
    
    #ifdef CRAYPAT
    PAT_record(PAT_STATE_OFF);
    #endif
    
    MPI_Finalize();
    return 0;
}
