//
//  main.cpp
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <iostream>
#include <omp.h>
#include <cassert>
#include <cstdlib>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

using namespace std;

//#define DEBUG

void tileTRD( const int MT, const int NT, TMatrix& A, TMatrix& T );

int main(int argc, const char * argv[])
{
	if (argc < 5)
	{
		cerr << "Usage: a.out [M] [N] [NB] [IB]\n";
		exit (1);
	}

	const int M =  atoi(argv[1]);  // n. of rows of the matrix
	const int N =  atoi(argv[2]);  // n. of columns of the matrix
	const int NB = atoi(argv[3]);  // tile size
	const int IB = atoi(argv[4]);  // inner blocking size

	assert( M >= N );
	assert( NB >= IB );

	#ifdef DEBUG
	cout << "M = " << M << ", N = " << N << ", NB = " << NB << ", IB = " << IB << endl;
	#endif

	//////////////////////////////////////////////////////////////////////
	// Definitions and Initialize
	TMatrix A(M,N,NB,NB,IB);

	const int MT = A.mt();
	const int NT = A.nt();

	// refered in workspace.c of PLASMA
	TMatrix T(MT*IB,NT*NB,IB,NB,IB);

	// Initialize matrix A
	A.Set_Rnd( 20170621 );

//	for (int i=0; i<MT; i++)
//		for (int j=0; j<NT; j++)
//			A(i,j)->Show_all();

	// Symmetrize
	for (int i=0; i<MT; i++)
		for (int j=i; j<NT; j++)
			for (int ii=0; ii<A(i,j)->m(); ii++)
				for (int jj=0; jj<A(i,j)->n(); jj++)
					A(j,i)->Set_Val(jj,ii,A(i,j)->Get_Val(ii,jj));

//	A.Show_all();
	// Definitions and Initialize　END
	//////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////

  // Timer start
  double time = omp_get_wtime();

  //////////////////////////////////////////////////////////////////////
  // tile QR variants
  tileTRD(MT,NT,A,T);
  //////////////////////////////////////////////////////////////////////

  // Timer stop
  time = omp_get_wtime() - time;
  cout << M << ", " << NB << ", " << IB << ", " << time << endl;

//	A.Show_all();
//	A(1,1)->Show_all();
//	A(2,2)->Show_all();
//	A(3,1)->Show_all();

  return EXIT_SUCCESS;
}
