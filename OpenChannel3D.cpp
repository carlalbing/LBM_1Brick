#include "OpenChannel3D.h"
#include "lattice_vars.h"
#include <fstream>
#include <string>
#include <sstream>

// for debugging
#include <iostream>

using namespace std;

OpenChannel3D::OpenChannel3D(const int rk, const int sz, const string input_file):
	rank(rk), size(sz)
{
	tag_d = 666; tag_u = 999;
	read_input_file(input_file);
	initialize_lattice_data();
	initialize_local_partition_variables();
	initialize_mpi_buffers();
	vtk_ts = 0;
}

OpenChannel3D::~OpenChannel3D(){
delete [] fEven;
delete [] fOdd;
delete [] inl;
delete [] onl;
delete [] snl;
delete [] u_bc;

delete [] ghost_in_m;
delete [] ghost_out_m;
delete [] ghost_in_p;
delete [] ghost_out_p;

delete [] rho_l;
delete [] ux_l;
delete [] uy_l;
delete [] uz_l;

}

void OpenChannel3D::write_data(MPI_Comm comm, bool isEven){

	string densityFileStub("density");
	string ux_FileStub("ux");
	string uy_FileStub("uy");
	string uz_FileStub("uz");
	stringstream ts_ind;
	string ts_ind_str;
	string fileSuffix(".b_dat");
	string density_fn, ux_fn, uy_fn, uz_fn;
	MPI_File fh_rho, fh_ux, fh_uy, fh_uz;
	MPI_Status mpi_s1, mpi_s2, mpi_s3,mpi_s4;

//	// create temporary data buffers
//	float * rho_l = new float[numEntries];
//	float * ux_l = new float[numEntries];
//	float * uy_l = new float[numEntries];
//	float * uz_l = new float[numEntries];

	float * fOut; 
	if (isEven){
		fOut = fEven;
	}else{
		fOut = fOdd;
	}

	#pragma omp parallel for collapse(3)
	for(int z = HALO;z<(totalSlices-HALO);z++){
	  for(int y = 0;y<Ny;y++){
	    for(int x = 0;x<Nx;x++){
	      int tid_l, tid_g;
	      float tmp_rho, tmp_ux, tmp_uy, tmp_uz;
	      tid_l = x+y*Nx+(z-HALO)*Nx*Ny;
	      tmp_rho = 0; tid_g = x+y*Nx+z*Nx*Ny;
	      tmp_ux = 0; tmp_uy = 0; tmp_uz = 0;
	      for(int spd=0;spd<numSpd;spd++){
		tmp_rho+=fOut[spd+tid_g*numSpd];
		tmp_ux+=ex[spd]*fOut[spd+tid_g*numSpd];
		tmp_uy+=ey[spd]*fOut[spd+tid_g*numSpd];
		tmp_uz+=ez[spd]*fOut[spd+tid_g*numSpd];
	      }
	      rho_l[tid_l]=tmp_rho;
	      if(snl[tid_g]==1){
		ux_l[tid_l]=0.; uy_l[tid_l]=0.; uz_l[tid_l]=0.;
	      }else{
		ux_l[tid_l]=tmp_ux*(1./tmp_rho);
		uy_l[tid_l]=tmp_uy*(1./tmp_rho);
		uz_l[tid_l]=tmp_uz*(1./tmp_rho);
	      }
	    }
	  }
	}

	// generate file names
	ts_ind << vtk_ts;
	density_fn = densityFileStub+ts_ind.str()+fileSuffix;
	ux_fn = ux_FileStub+ts_ind.str()+fileSuffix;
	uy_fn = uy_FileStub+ts_ind.str()+fileSuffix;
	uz_fn = uz_FileStub+ts_ind.str()+fileSuffix;

	// open MPI file for parallel IO
	MPI_File_open(comm,(char*)density_fn.c_str(),
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_rho);

	MPI_File_open(comm,(char*)ux_fn.c_str(),
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_ux);

	MPI_File_open(comm,(char*)uy_fn.c_str(),
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uy);

	MPI_File_open(comm,(char*)uz_fn.c_str(),
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh_uz);

	//write your chunk of data
	MPI_File_write_at(fh_rho,offset,rho_l,numEntries,MPI_FLOAT,&mpi_s1);
	MPI_File_write_at(fh_ux,offset,ux_l,numEntries,MPI_FLOAT,&mpi_s2);
	MPI_File_write_at(fh_uy,offset,uy_l,numEntries,MPI_FLOAT,&mpi_s3);
	MPI_File_write_at(fh_uz,offset,uz_l,numEntries,MPI_FLOAT,&mpi_s4);

	//close the files
	MPI_File_close(&fh_rho);
	MPI_File_close(&fh_ux);
	MPI_File_close(&fh_uy);
	MPI_File_close(&fh_uz);

	vtk_ts++; // increment the dump counter...


}

void OpenChannel3D::D3Q15_process_slices(bool isEven, const int firstSlice,
		  const int lastSlice){

  // this monstrosity needs to be change into something more simple and clear.
  // it performs on the GPU but is rather unmaintainable.

  const float * fIn;
  float * fOut;

  if(isEven){
    fIn = fEven; fOut = fOdd;
  }else{
    fIn = fOdd; fOut = fEven;
  }


  //Nz=lastSlice-firstSlice;
  const int numSpd=15;
  #pragma omp parallel for collapse(3)
  for(int Z=firstSlice;Z<lastSlice;Z++){
    for(int Y=0;Y<Ny;Y++){
      for(int X=0;X<Nx;X++){
	float f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14;
	float cu,rho,ux,uy,uz,fEq,dz;
	int X_t,Y_t,Z_t,tid_t,tid;

	tid=X+Y*Nx+Z*Nx*Ny;

	//load the data into registers
	// f0=fIn[tid]; f1=fIn[Nx*Ny*Nz+tid];
	// f2=fIn[2*Nx*Ny*Nz+tid]; f3=fIn[3*Nx*Ny*Nz+tid];
	// f4=fIn[4*Nx*Ny*Nz+tid]; f5=fIn[5*Nx*Ny*Nz+tid];
	// f6=fIn[6*Nx*Ny*Nz+tid]; f7=fIn[7*Nx*Ny*Nz+tid];
	// f8=fIn[8*Nx*Ny*Nz+tid]; f9=fIn[9*Nx*Ny*Nz+tid];
	// f10=fIn[10*Nx*Ny*Nz+tid]; f11=fIn[11*Nx*Ny*Nz+tid];
	// f12=fIn[12*Nx*Ny*Nz+tid]; f13=fIn[13*Nx*Ny*Nz+tid];
	// f14=fIn[14*Nx*Ny*Nz+tid];
	f0=fIn[tid*numSpd]; f1=fIn[tid*numSpd+1];
	f2=fIn[tid*numSpd+2]; f3=fIn[tid*numSpd+3];
	f4=fIn[tid*numSpd+4]; f5=fIn[tid*numSpd+5];
	f6=fIn[tid*numSpd+6]; f7=fIn[tid*numSpd+7];
	f8=fIn[tid*numSpd+8]; f9=fIn[tid*numSpd+9];
	f10=fIn[tid*numSpd+10]; f11=fIn[tid*numSpd+11];
	f12=fIn[tid*numSpd+12]; f13=fIn[tid*numSpd+13];
	f14=fIn[tid*numSpd+14];

	//compute density
	rho = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14;
	ux=f1-f2+f7-f8+f9-f10+f11-f12+f13-f14; ux/=rho;
	uy=f3-f4+f7+f8-f9-f10+f11+f12-f13-f14; uy/=rho;
	uz=f5-f6+f7+f8+f9+f10-f11-f12-f13-f14; uz/=rho;
	
	//if it's on the inl or onl, update

	if((inl[tid]==1)||(onl[tid]==1)){

	  dz=u_bc[tid]-uz;
	  //speed 1 ex=1 ey=ez=0. w=1./9.
	  cu=3.*(1.)*(-ux);
	  f1+=(1./9.)*rho*cu;

	  //speed 2 ex=-1 ey=ez=0. w=1./9.
	  cu=3.*(-1.)*(-ux);
	  f2+=(1./9.)*rho*cu;

	  //speed 3 ey=1; ex=ez=0; w=1./9.
	  cu=3.*(1.)*(-uy);
	  f3+=(1./9.)*rho*cu;

	  //speed 4 ey=-1; ex=ez=0; w=1./9.
	  cu=3.*(-1.)*(-uy);
	  f4+=(1./9.)*rho*cu;

	  //speed 5 ex=ey=0; ez=1; w=1./9.
	  cu=3.*(1.)*(dz);
	  f5+=(1./9.)*rho*cu;

	  //speed 6 ex=ey=0; ez=-1; w=1./9.
	  cu=3.*(-1.)*(dz);
	  f6+=(1./9.)*rho*cu;

	  //speed 7 ex=ey=ez=1; w=1./72.
	  cu=3.*((1.)*-ux+(1.)*(-uy)+(1.)*dz);
	  f7+=(1./72.)*rho*cu;

	  //speed 8 ex=-1 ey=ez=1; w=1./72.
	  cu=3.*((-1.)*-ux+(1.)*(-uy)+(1.)*dz);
	  f8+=(1./72.)*rho*cu;

	  //speed 9 ex=1 ey=-1 ez=1
	  cu=3.0*((1.)*-ux+(-1.)*(-uy)+(1.)*dz);
	  f9+=(1./72.)*rho*cu;

	  //speed 10 ex=-1 ey=-1 ez=1
	  cu=3.0*((-1.)*-ux+(-1.)*(-uy)+(1.)*dz);
	  f10+=(1./72.)*rho*cu;

	  //speed 11 ex=1 ey=1 ez=-1
	  cu=3.0*((1.)*-ux +(1.)*(-uy)+(-1.)*dz);
	  f11+=(1./72.)*rho*cu;

	  //speed 12 ex=-1 ey=1 ez=-1
	  cu=3.0*((-1.)*-ux+(1.)*(-uy)+(-1.)*dz);
	  f12+=(1./72.)*rho*cu;

	  //speed 13 ex=1 ey=-1 ez=-1 w=1./72.
	  cu=3.0*((1.)*-ux+(-1.)*(-uy)+(-1.)*dz);
	  f13+=(1./72.)*rho*cu;
      
	  //speed 14 ex=ey=ez=-1 w=1./72.
	  cu=3.0*((-1.)*-ux + (-1.)*(-uy) +(-1.)*dz);
	  f14+=(1./72.)*rho*cu;

	  ux=0.; uy=0.; uz=u_bc[tid];
	}

	if(snl[tid]==1){
	  // 1--2
	  cu=f2; f2=f1; f1=cu;
	  //3--4
	  cu=f4; f4=f3; f3=cu;
	  //5--6
	  cu=f6; f6=f5; f5=cu;
	  //7--14
	  cu=f14; f14=f7; f7=cu;
	  //8--13
	  cu=f13; f13=f8; f8=cu;
	  //9--12
	  cu=f12; f12=f9; f9=cu;
	  //10--11
	  cu=f11; f11=f10; f10=cu;


	}else{
	  fEq=rho*(2./9.)*(1.-1.5*(ux*ux+uy*uy+uz*uz));
	  f0=f0-omega*(f0-fEq);

	  //speed 1 ex=1 ey=ez=0 w=1./9.
	  cu=3.*(1.*ux);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f1=f1-omega*(f1-fEq);

	  //speed 2 ex=-1 ey=ez=0 w=1./9.
	  cu=3.*((-1.)*ux);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f2=f2-omega*(f2-fEq);

	  //speed 3 ex=0 ey=1 ez=0 w=1./9.
	  cu=3.*(1.*uy);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f3=f3-omega*(f3-fEq);

	  //speed 4 ex=0 ey=-1 ez=0 w=1./9.
	  cu=3.*(-1.*uy);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f4=f4-omega*(f4-fEq);

	  //speed 5 ex=ey=0 ez=1 w=1./9.
	  cu=3.*(1.*uz);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f5=f5-omega*(f5-fEq);

	  //speed 6 ex=ey=0 ez=-1 w=1./9.
	  cu=3.*(-1.*uz);
	  fEq=rho*(1./9.)*(1.+cu+0.5*(cu*cu)-
			   1.5*(ux*ux+uy*uy+uz*uz));
	  f6=f6-omega*(f6-fEq);

	  //speed 7 ex=ey=ez=1 w=1./72.
	  cu=3.*(ux+uy+uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f7=f7-omega*(f7-fEq);

	  //speed 8 ex=-1 ey=ez=1 w=1./72.
	  cu=3.*(-ux+uy+uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f8=f8-omega*(f8-fEq);

	  //speed 9 ex=1 ey=-1 ez=1 w=1./72.
	  cu=3.*(ux-uy+uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f9=f9-omega*(f9-fEq);

	  //speed 10 ex=-1 ey=-1 ez=1 w=1/72
	  cu=3.*(-ux-uy+uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f10=f10-omega*(f10-fEq);

	  //speed 11 ex=1 ey=1 ez=-1 w=1/72
	  cu=3.*(ux+uy-uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f11=f11-omega*(f11-fEq);

	  //speed 12 ex=-1 ey=1 ez=-1 w=1/72
	  cu=3.*(-ux+uy-uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f12=f12-omega*(f12-fEq);

	  //speed 13 ex=1 ey=ez=-1 w=1/72
	  cu=3.*(ux-uy-uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f13=f13-omega*(f13-fEq);

	  //speed 14 ex=ey=ez=-1 w=1/72
	  cu=3.*(-ux-uy-uz);
	  fEq=rho*(1./72.)*(1.+cu+0.5*(cu*cu)-
			    1.5*(ux*ux+uy*uy+uz*uz));
	  f14=f14-omega*(f14-fEq);



	}

	//speed 0 ex=ey=ez=0
	//fOut[tid]=f0;
	fOut[tid*numSpd]=f0;

	//speed 1 ex=1 ey=ez=0
	X_t=X+1; Y_t=Y; Z_t=Z;
	if(X_t==Nx) X_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[Nx*Ny*Nz+tid_t]=f1;
	fOut[tid_t*numSpd+1]=f1;

	//speed 2 ex=-1 ey=ez=0;
	X_t=X-1; Y_t=Y; Z_t=Z;
	if(X_t<0) X_t=(Nx-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[2*Nx*Ny*Nz+tid_t]=f2;
	fOut[tid_t*numSpd+2]=f2;

	//speed 3 ex=0 ey=1 ez=0
	X_t=X; Y_t=Y+1; Z_t=Z;
	if(Y_t==Ny) Y_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[3*Nx*Ny*Nz+tid_t]=f3;
	fOut[tid_t*numSpd+3]=f3;

	//speed 4 ex=0 ey=-1 ez=0
	X_t=X; Y_t=Y-1; Z_t=Z;
	if(Y_t<0) Y_t=(Ny-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	///	fOut[4*Nx*Ny*Nz+tid_t]=f4;
	fOut[tid_t*numSpd+4]=f4;


	//speed 5 ex=ey=0 ez=1
	X_t=X; Y_t=Y; Z_t=Z+1;
	//	if(Z_t==Nz) Z_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//fOut[5*Nx*Ny*Nz+tid_t]=f5;
	fOut[tid_t*numSpd+5]=f5;

	//speed 6 ex=ey=0 ez=-1
	X_t=X; Y_t=Y; Z_t=Z-1;
	//	if(Z_t<0) Z_t=(Nz-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[6*Nx*Ny*Nz+tid_t]=f6;
	fOut[tid_t*numSpd+6]=f6;

	//speed 7 ex=ey=ez=1
	X_t=X+1; Y_t=Y+1; Z_t=Z+1;
	if(X_t==Nx) X_t=0;
	if(Y_t==Ny) Y_t=0;
	//	if(Z_t==Nz) Z_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[7*Nx*Ny*Nz+tid_t]=f7;
	fOut[tid_t*numSpd+7]=f7;

	//speed 8 ex=-1 ey=1 ez=1
	X_t=X-1; Y_t=Y+1; Z_t=Z+1;
	if(X_t<0) X_t=(Nx-1);
	if(Y_t==Ny) Y_t=0;
	//	if(Z_t==Nz) Z_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[8*Nx*Ny*Nz+tid_t]=f8;
	fOut[tid_t*numSpd+8]=f8;

	//speed 9 ex=1 ey=-1 ez=1
	X_t=X+1; Y_t=Y-1; Z_t=Z+1;
	if(X_t==Nx) X_t=0;
	if(Y_t<0) Y_t=(Ny-1);
	//	if(Z_t==Nz) Z_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[9*Nx*Ny*Nz+tid_t]=f9;
	fOut[tid_t*numSpd+9]=f9;

	//speed 10 ex=-1 ey=-1 ez=1
	X_t=X-1; Y_t=Y-1; Z_t=Z+1;
	if(X_t<0) X_t=(Nx-1);
	if(Y_t<0) Y_t=(Ny-1);
	//	if(Z_t==Nz) Z_t=0;
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[10*Nx*Ny*Nz+tid_t]=f10;
	fOut[tid_t*numSpd+10]=f10;

	//speed 11 ex=1 ey=1 ez=-1
	X_t=X+1; Y_t=Y+1; Z_t=Z-1;
	if(X_t==Nx) X_t=0;
	if(Y_t==Ny) Y_t=0;
	//	if(Z_t<0) Z_t=(Nz-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[11*Nx*Ny*Nz+tid_t]=f11;
	fOut[tid_t*numSpd+11]=f11;

	//speed 12 ex=-1 ey=1 ez=-1
	X_t=X-1; Y_t=Y+1; Z_t=Z-1;
	if(X_t<0) X_t=(Nx-1);
	if(Y_t==Ny) Y_t=0;
	//	if(Z_t<0) Z_t=(Nz-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[12*Nx*Ny*Nz+tid_t]=f12;
	fOut[tid_t*numSpd+12]=f12;

	//speed 13 ex=1 ey=-1 ez=-1
	X_t=X+1; Y_t=Y-1; Z_t=Z-1;
	if(X_t==Nx) X_t=0;
	if(Y_t<0) Y_t=(Ny-1);
	//	if(Z_t<0) Z_t=(Nz-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[13*Nx*Ny*Nz+tid_t]=f13;
	fOut[tid_t*numSpd+13]=f13;

	//speed 14 ex=ey=ez=-1
	X_t=X-1; Y_t=Y-1; Z_t=Z-1;
	if(X_t<0) X_t=(Nx-1);
	if(Y_t<0) Y_t=(Ny-1);
	//	if(Z_t<0) Z_t=(Nz-1);
	tid_t=X_t+Y_t*Nx+Z_t*Nx*Ny;
	//	fOut[14*Nx*Ny*Nz+tid_t]=f14;

	fOut[tid_t*numSpd+14]=f14;
      }
    }
  }


}

void OpenChannel3D::stream_out_collect(const float * fIn_b, float * buff_out,
		const int numStreamSpeeds,
		const int * streamSpeeds){
	int numSpeeds = numSpd;
	int tid_l, stream_spd;
	for(int z=0;z<HALO;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				for(int spd=0;spd<numStreamSpeeds;spd++){
					tid_l = x+y*Nx+z*Nx*Ny;
					stream_spd=streamSpeeds[spd];
					buff_out[tid_l*numStreamSpeeds+spd]=fIn_b[tid_l*numSpeeds+stream_spd];
				}
			}
		}
	}
}

void OpenChannel3D::stream_in_distribute(float * fIn_b,
		const float * buff_in, const int numStreamSpeeds,
		const int * streamSpeeds){

	int tid_l, stream_spd;
	for(int z=0;z<HALO;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				for(int spd=0;spd<numStreamSpeeds;spd++){
					tid_l=x+y*Nx+z*Nx*Ny;
					stream_spd=streamSpeeds[spd];
					fIn_b[tid_l*numSpd+stream_spd]=buff_in[tid_l*numStreamSpeeds+spd];
				}

			}
		}
	}



}
void OpenChannel3D::take_lbm_timestep(bool isEven, MPI_Comm comm){


	// collide and stream lower boundary slices
	D3Q15_process_slices(isEven,HALO,HALO+1);

	// collect data in ghost_m_out
	float * fIn_b;
	if(isEven){
		fIn_b = ghost_out_odd_m;

	}else{
		fIn_b = ghost_out_even_m;
	}

	// cout << "rank: " << rank << "in take_lbm_timestep: " << endl;
	// cout << "rank: " << rank << "fIn_b (ghost_out_[odd/even]_m buffer) = " << fIn_b << endl;

	stream_out_collect(fIn_b,ghost_out_m,numMspeeds,Mspeeds);
	// begin communication to ghost_p_in
	MPI_Isend(ghost_out_m,numHALO,MPI_FLOAT,nd_m,tag_d,comm, &rq_out1);
	MPI_Irecv(ghost_in_p,numHALO,MPI_FLOAT,nd_p,tag_d,comm,&rq_in1);
	// collide and stream upper boundary slices
	D3Q15_process_slices(isEven,totalSlices-2*HALO,totalSlices-HALO);
	// collect data in ghost_p_out
	if(isEven){
		fIn_b = ghost_out_odd_p;

	}else{
		fIn_b = ghost_out_even_p;
	}

	//	cout << "rank: " << rank << "(ghost_out_[odd/even]_p) = " << fIn_b << endl;

	stream_out_collect(fIn_b,ghost_out_p,numPspeeds,Pspeeds);
	// begin communication to ghost_m_in
	MPI_Isend(ghost_out_p,numHALO,MPI_FLOAT,nd_p,tag_u,comm,&rq_out2);
	MPI_Irecv(ghost_in_m,numHALO,MPI_FLOAT,nd_m,tag_u,comm,&rq_in2);
	// collide and stream interior lattice points
	D3Q15_process_slices(isEven,HALO+1,totalSlices-2*HALO);
	// ensure communication of boundary lattice points is complete
	MPI_Wait(&rq_in1,&stat);
	MPI_Wait(&rq_in2,&stat);
	// copy data from ghost_in_p to Mspeeds on P boundary points
	if(isEven){
		fIn_b = ghost_in_odd_p;
	}else{
		fIn_b = ghost_in_even_p;
	}
	stream_in_distribute(fIn_b,ghost_in_p,numMspeeds,Mspeeds);

	//	cout << "rank: " << rank << "ghost_in_[odd/even]_p = " << fIn_b << endl;

	// copy data from ghost_m_in to Pspeeds on M boundary points
	if(isEven){
		fIn_b = ghost_in_odd_m;
	}else{
		fIn_b = ghost_in_even_m;
	}
	stream_in_distribute(fIn_b,ghost_in_m,numPspeeds,Pspeeds);
	//cout << "rank: " << rank << "ghost_in_[odd/even]_m = " << fIn_b << endl;

}

void OpenChannel3D::initialize_mpi_buffers(){
	// integer for arithmetic into pointer arrays for nodes I need to
	// communicate to neighbors.
	firstNdm = Nx*Ny*HALO;
	lastNdm = Nx*Ny*(HALO+1);
	firstNdp = nnodes-2*(Nx*Ny*HALO);
	lastNdp = nnodes-(Nx*Ny*HALO);

	// actual pointers
	ghost_out_odd_p = fOdd +(Nx*Ny*(totalSlices-HALO)*numSpd);
	ghost_out_odd_m = fOdd;
	ghost_in_odd_p = fOdd+(Nx*Ny*numSpd*numMySlices);
	ghost_in_odd_m = fOdd+(Nx*Ny*numSpd*HALO);

	// cout << "rank: " << rank << "In initialize MPI buffers: " << endl;
	// cout << "rank: " << rank << "ghost_out_odd_p = " << ghost_out_odd_p << endl;
	// cout << "rank: " << rank << "ghost_out_odd_m = " << ghost_out_odd_m << endl;
	// cout << "rank: " << rank << "ghost_in_odd_p = " << ghost_in_odd_p << endl;
	// cout << "rank: " << rank << "ghost_in_odd_m = " << ghost_in_odd_m << endl;

	ghost_out_even_p = fEven+(Nx*Ny*(totalSlices-HALO)*numSpd);
	ghost_out_even_m = fEven;
	ghost_in_even_p = fEven+(Nx*Ny*numSpd*numMySlices);
	ghost_in_even_m = fEven+(Nx*Ny*numSpd*HALO);

	// cout << "rank: " << rank << "ghost_out_even_p = " << ghost_out_even_p << endl;
	// cout << "rank: " << rank << "ghost_out_even_m = " << ghost_out_even_m << endl;
	// cout << "rank: " << rank << "ghost_in_even_p = " << ghost_in_even_p << endl;
	// cout << "rank: " << rank << "ghost_in_even_m = " << ghost_in_even_m << endl;


	numHALO = (Nx*Ny*numPspeeds*HALO);

	ghost_in_m = new float[Nx*Ny*numMspeeds*HALO];
	ghost_out_m = new float[Nx*Ny*numMspeeds*HALO];
	ghost_in_p = new float[Nx*Ny*numPspeeds*HALO];
	ghost_out_p = new float[Nx*Ny*numPspeeds*HALO];
	offset = firstSlice*Nx*Ny*sizeof(float);
	numEntries = numMySlices*Nx*Ny;

	rho_l = new float[numEntries];
	ux_l = new float[numEntries];
	uy_l = new float[numEntries];
	uz_l = new float[numEntries];



}

void OpenChannel3D::initialize_local_partition_variables(){
	numMySlices = Nz/size;
	if(rank<(Nz%size))
		numMySlices+=1;

	firstSlice=(Nz/size)*rank;
	if((Nz%size)<rank){
		firstSlice+=(Nz%size);
	}else{
		firstSlice+=rank;
	}
	lastSlice = firstSlice+numMySlices-1;
	totalSlices=numMySlices+2*HALO;
	nnodes = totalSlices*Nx*Ny;

	fEven = new float[nnodes*numSpd];
	fOdd = new float[nnodes*numSpd];
	snl = new int[nnodes];
	inl = new int[nnodes];
	onl = new int[nnodes];
	u_bc = new float[nnodes];


	// cout << "rank: " << rank << "In initialize_local_partition_variables:" << endl << " fEven = " << fEven << endl;
	// cout << "rank: " << rank << " fOdd = " << fOdd << endl;


	nd_m = rank-1;
	if(nd_m<0)
		nd_m=(size-1);

	nd_p = rank+1;
	if(nd_p==size)
		nd_p=0;


	// initialize the variables
	int tid;
	for(int z=0; z<totalSlices;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				tid=x+y*Nx+z*Nx*Ny;
				for(int spd=0;spd<numSpd;spd++){
					fEven[spd+tid*numSpd]=rho_lbm*w[spd];
					fOdd[spd+tid*numSpd]=rho_lbm*w[spd];
				}
			}
		}
	}

	// initialize all node lists to zero...
	int tid_l, x,z, y;
	for(int nd=0;nd<nnodes;nd++){
		inl[nd]=0;
		onl[nd]=0;
		snl[nd]=0;
	}
	// initialize the solid node list
	y=0;
	for(int z=0;z<totalSlices;z++){
		for(int x=0;x<Nx;x++){
			tid_l=x+y*Nx+z*Nx*Ny;
			snl[tid_l]=1;
		}
	}
	y=(Ny-1);
	for(int z=0;z<totalSlices;z++){
		for(int x=0;x<Nx;x++){
			tid_l=x+y*Nx+z*Nx*Ny;
			snl[tid_l]=1;
		}
	}

	// near and far wall
	for(int z=0;z<totalSlices;z++){
		for(int y=0;y<Ny;y++){
			x = 0;
			tid_l = x+y*Nx+z*Nx*Ny;
			snl[tid_l]=1;
			x = (Nx-1);
			tid_l = x+y*Nx+z*Nx*Ny;
			snl[tid_l]=1;
		}
	}



	// Due to how the domain is partitioned, the inlet nodes
	// are all assigned to rank 0 process
	if(rank==0){
		z=1; // to account for the HALO nodes on rank 0, this is z=1, not z=0
		for(int y=1;y<(Ny-1);y++){//<-- skip top and bottom
			for(int x=1;x<(Nx-1);x++){
				tid_l = x+y*Nx+z*Nx*Ny;
				inl[tid_l]=1;
			}
		}
	}
	// again, due to the partitioning strategy, all of the outlet nodes
	// are on the (size-1) partition
	if(rank==(size-1)){
		z=totalSlices-1; // again, to account for the HALO at the outlet, this is z = totalSlices-1, not z=totalSlices...
		for(int y=1;y<(Ny-1);y++){
			for(int x=1;x<(Nx-1);x++){
				tid_l=x+y*Nx+z*Nx*Ny;
				inl[tid_l]=1;
			}
		}
	}


	// rank 0 partition has the outlet slice on the left halo (due to logical periodicity)
	if(rank==0){
		z=0; // inlet halo is he outlet slice, so z=0 on the rank 0 process is the outlet slice.
		for(int y=1;y<(Ny-1);y++){
			for(int x=1;x<(Nx-1);x++){
				tid_l=x+y*Nx+z*Nx*Ny;
				onl[tid_l]=1;
			}
		}
	}
	// rank (size-1) has the inlet slice on the right halo. (due to logical periodicity)
	if(rank==(size-1)){
		z=totalSlices-2; // because of the HALO, the outlet nodes are at totalSlices-2, not totalSlices-1
		for(int y=1;y<(Ny-1);y++){
			for(int x=1;x<(Nx-1);x++){
				tid_l=x+y*Nx+z*Nx*Ny;
				onl[tid_l]=1;
			}
		}

	}

	// initialize u_bc
	float b = ((float)Ny-1.)/2.;
	float h;
	for(int z=0;z<totalSlices;z++){
		for(int y=0;y<Ny;y++){
			for(int x=0;x<Nx;x++){
				tid=x+y*Nx+z*Nx*Ny;
				if((inl[tid]==1)|(onl[tid]==1)){
					h=((float)y-b)/b;
					u_bc[tid]=umax_lbm*(1.-(h*h));
				}else{
					u_bc[tid]=0.;
				}
			}
		}
	}

	// initialize outlet density
	rho_out = rho_lbm;


}

void OpenChannel3D::read_input_file(const string input_file){
	ifstream input_params(input_file.c_str(),ios::in);
	input_params >> LatticeType;
	input_params >> Num_ts;
	input_params >> ts_rep_freq;
	input_params >> plot_freq;
	/*input_params >> obst_type;
	input_params >> obst_param1;
	input_params >> obst_param2;
	input_params >> obst_param3;
	input_params >> obst_param4;*/
	input_params >> rho_lbm;
	input_params >> umax_lbm;
	input_params >> omega;
	input_params >> Nx;
	input_params >> Ny;
	input_params >> Nz;
	input_params.close();


}

void OpenChannel3D::initialize_lattice_data(){

	switch(LatticeType){
	case(1):
				numSpd=15;
	ex = ex15;
	ey = ey15;
	ez = ez15;
	w = w15;
	numPspeeds = numPspeedsD3Q15;
	numMspeeds = numMspeedsD3Q15;
	Mspeeds = MspeedsD3Q15;
	Pspeeds = PspeedsD3Q15;


	}


}
