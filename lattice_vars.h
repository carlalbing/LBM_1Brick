#ifndef LATTICE_VARS_H
#define LATTICE_VARS_H

//------------------- D2Q9 ------------------------------
 float ex9[9]={0,1,0,-1,0,1,-1,-1,1};
 float ey9[9]={0,0,1,0,-1,1,1,-1,-1};
 float w9[9]={4.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f};
 int bb9[9]={0,3,4,1,2,7,8,5,6};

 float M9[9*9]={1,1,1,1,1,1,1,1,1,
	      -4,-1,-1,-1,-1,2,2,2,2,
	      4,-2,-2,-2,-2,1,1,1,1,
	      0,1,0,-1,0,1,-1,-1,1,
	      0,-2,0,2,0,1,-1,-1,1,
	      0,0,1,0,-1,1,1,-1,-1,
	      0,0,-2,0,2,1,1,-1,-1,
	      0,1,-1,1,-1,0,0,0,0,
	      0,0,0,0,0,1,-1,1,-1};

//-------- D3Q15 ---------------------------------------
 float ex15[15]={0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1};
 float ey15[15]={0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1};
 float ez15[15]={0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1};
 float w15[15]={2.f/9.f,1.f/9.f,1.f/9,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,
	       1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f,
	       1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f};
 int bb15[15]={0,2,1,4,3,6,5,14,13,12,11,10,9,8,7};
 float M15[15*15]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			-2,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,
			16,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,
			0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,
			0,-4,4,0,0,0,0,1,-1,1,-1,1,-1,1,-1,
			0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1,
			0,0,0,-4,4,0,0,1,1,-1,-1,1,1,-1,-1,
			0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1,
			0,0,0,0,0,-4,4,1,1,1,1,-1,-1,-1,-1,
			0,2,2,-1,-1,-1,-1,0,0,0,0,0,0,0,0,
			0,0,0,1,1,-1,-1,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,1,-1,-1,1,1,-1,-1,1,
			0,0,0,0,0,0,0,1,1,-1,-1,-1,-1,1,1,
			0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,
			0,0,0,0,0,0,0,1,-1,-1,1,-1,1,1,-1};

int numPspeedsD3Q15=5;
int PspeedsD3Q15[5]={5,7,8,9,10};
int numMspeedsD3Q15=5;
int MspeedsD3Q15[5]={6,11,12,13,14};
//------------------------------------------------------

//---------- D3Q19 ------------------------------------
 float ex19[19]={0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0};
 float ey19[19]={0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1};
 float ez19[19]={0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1};
 float w19[19]={3.f/9.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,
	       1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,
	       1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f};
 int bb19[19]={0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
 float M19[19*19]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			-30, -11, -11, -11, -11, -11, -11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
			12, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1,0, 0, 0, 0,
			0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0,
			0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, -4, 4, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1,
			0, 2, 2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2,
			0, -4, -4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2,
			0, 0, 0, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
			0, 0, 0, -2, -2, 2, 2, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0,0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, -1, -1, 1, 1, 0, 0, 0, 0, 1, -1, 1, -1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1};

//----------------------------------------------------

//--------- D3Q27 ----------------------------------------
 float ex27[27]={0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1,0,0,1,1,1,1};
 float ey27[27]={0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,1,-1,0,0,1,1,1,1,-1,-1};
 float ez27[27]={0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1,1,-1,1,-1};
 float w27[27]={8.f/27.f,2.f/27.f,2.f/27.f,2.f/27.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,
	       1.f/54.f,1.f/216.f,1.f/216,1.f/216.f,1.f/216.f,2.f/27.f,2.f/27.f,
	       2.f/27.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,
		1.f/54.f,1.f/216.f,1.f/216,1.f/216.f,1.f/216.f};
 int bb27[27]={0,14,15,16,17,18,19,20,21,22,23,24,25,26,
	      1,2,3,4,5,6,7,8,9,10,11,12,13};


//--------------------------------------------------------


#endif
