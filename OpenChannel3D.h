#ifndef OPENCHANNEL3D_H
#define OPENCHANNEL3D_H

#include <mpi.h>
#include <string>
#define HALO 1

using namespace std;

class OpenChannel3D{
public:

	// MPI data
	const int rank,size;
	MPI_Status stat;
	MPI_Request rq_in1,rq_in2;
	MPI_Request rq_out1,rq_out2;
	int tag_d, tag_u;

	//input data
	int LatticeType;
	int Num_ts;
	int ts_rep_freq;
	int plot_freq;
	/*int obst_type;
	float obst_param1;
	float obst_param2;
	float obst_param3;
	float obst_param4;*/
	float rho_lbm;
	float umax_lbm;
	float omega;
	int Nx;
	int Ny;
	int Nz;

	// lattice data
	int numSpd;
	float * ex;
	float * ey;
	float * ez;
	float * w;
	int numPspeeds;
	int * Pspeeds;
	int numMspeeds;
	int * Mspeeds;

	// local partition variables
	int numMySlices,numMyNodes, firstSlice, lastSlice, totalSlices, nnodes;
	float * fEven, *fOdd, * u_bc;
	float rho_out;
	int nd_m, nd_p;
	int * snl, * inl, * onl;

	// MPI communication info
	int firstNdm, lastNdm, firstNdp, lastNdp, numHALO;
	// pointers to boundary slice info
	float * ghost_out_odd_p, * ghost_out_odd_m, * ghost_in_odd_p, * ghost_in_odd_m;
	float * ghost_out_even_p, * ghost_out_even_m, * ghost_in_even_p, * ghost_in_even_m;
	// data buffers
	float * ghost_in_m, * ghost_out_m, * ghost_in_p, * ghost_out_p;

	// visualization dump count variable
	int vtk_ts;
	// visualization file offset
	int offset, numEntries;
	float * rho_l, *ux_l, *uy_l, *uz_l;

	OpenChannel3D(const int rank, const int size, const string input_file);
  ~OpenChannel3D();
  void write_data(MPI_Comm comm, bool isEven);
  void take_lbm_timestep(bool isEven, MPI_Comm comm);

 private:
  void read_input_file(const string input_file);
  void initialize_lattice_data();
  void initialize_local_partition_variables();
  void initialize_mpi_buffers();
  void D3Q15_process_slices(bool isEven, const int firstSlice,
		  const int lastSlice);
  void stream_out_collect(const float * fIn_b, float * buff_out,
			const int numStreamSpeeds,
			const int * streamSpeeds);
  void stream_in_distribute(float * fIn_b,
			const float * buff_in, const int numStreamSpeeds,
			const int * streamSpeeds);






};

#endif
