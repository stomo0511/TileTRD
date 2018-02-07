//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

//#define COUT
//#define ANIM

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <omp.h>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>
#include <BMatrix.hpp>

using namespace std;

void tileCopy( BMatrix* A, BMatrix& B )
{
	assert( A->m() == B.m() );
	assert( A->m() == B.n() );
	assert( A->ib() == B.ib() );

	for (int i=0; i<A->m(); i++)
		for (int j=0; j<A->n(); j++)
			B.Set_Val(i,j,A->Get_Val(i,j));
}

void tileTransCopy( BMatrix* A, BMatrix& B )
{
	assert( A->m() == B.m() );
	assert( A->m() == B.n() );
	assert( A->ib() == B.ib() );

	for (int i=0; i<A->m(); i++)
		for (int j=0; j<A->n(); j++)
			B.Set_Val(j,i,A->Get_Val(i,j));
}

void tileTransCopy( BMatrix A, BMatrix* B )
{
	assert( A.m() == B->m() );
	assert( A.m() == B->n() );
	assert( A.ib() == B->ib() );

	for (int i=0; i<A.m(); i++)
		for (int j=0; j<A.n(); j++)
			B->Set_Val(j,i,A.Get_Val(i,j));
}

void tileTRD( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	// Progress table
	int **Ap, **Tp;

	Ap = (int **)malloc( sizeof(int*) * MT);
	for (int i=0; i<MT; i++)
		Ap[i] = (int *)malloc( sizeof(int) * NT);

	Tp = (int **)malloc( sizeof(int*) * MT);
	for (int i=0; i<MT; i++)
		Tp[i] = (int *)malloc( sizeof(int) * NT);

	const int NB = A(0,0)->m();
	const int IB = A(0,0)->ib();
	BMatrix tmp(NB,NB,IB);

	#ifdef ANIM
	cout << "Kernel,Ii,Ij,Ik,Time\n";
	#endif

	double ttime = omp_get_wtime();

	#pragma omp parallel firstprivate(ttime)
	{
		#pragma omp single
		{
			for (int tk=0; tk < min(MT,NT)-1; tk++ )
			{
				// (a)
				#pragma omp task depend (inout:Ap[tk+1][tk]) depend (out:Tp[tk+1][tk])
				{
					GEQRT( A(tk+1,tk), T(tk+1,tk) );
					#ifdef COUT
					#pragma omp critical
					cout << "GEQRT(" << tk+1 << "," << tk << "," << tk << ")\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "GL," << tk+1 << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

				// (b)
				#pragma omp task depend(in:Ap[tk+1][tk], Tp[tk+1][tk]) depend(inout:Ap[tk+1][tk+1])
				{
					LARFB( PlasmaLeft, PlasmaTrans,   A(tk+1,tk), T(tk+1,tk), A(tk+1,tk+1) );
					#ifdef COUT
					#pragma omp critical
					cout << "LARFB_L(" << tk+1 << "," << tk+1 << "," << tk << ")\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "LL," << tk+1 << "," << tk+1 << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

				// (c)
				for (int ti=tk+1; ti < MT; ti++)
				{
					#pragma omp task depend(in:Ap[tk+1][tk], Tp[tk+1][tk]) depend(inout:Ap[ti][tk+1])
					{
						LARFB( PlasmaRight, PlasmaNoTrans, A(tk+1,tk), T(tk+1,tk), A(ti,tk+1) );
						#ifdef COUT
						#pragma omp critical
						cout << "LARFB_R(" << ti << "," << tk+1 << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "LR," << ti << "," << tk+1 << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}
				} // i-LOOP END

				for (int ti=tk+2; ti < MT; ti++)
				{
					// (d)
					#pragma omp task depend(inout:Ap[tk+1][tk], Ap[ti][tk]) depend(out:Tp[ti][tk])
					{
						TSQRT( A(tk+1,tk), A(ti,tk), T(ti,tk) );
						#ifdef COUT
						#pragma omp critical
						cout << "TSQRT( (" << tk+1 << "," << tk << "," << tk << "), (" << ti << "," << tk << "," << tk << ") )\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "TL," << ti << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					//////////////////////////////////
					// (e) tj=tk+1
					#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) depend(inout: Ap[tk+1][tk+1], Ap[ti][tk+1]) depend(out: tmp)
					{
						tileTransCopy(A(ti,tk+1),tmp);  // tmp <- A(ti,tk+1)^T
						#ifdef COUT
						#pragma omp critical
						cout << "tmp <- (" << ti << "," << tk+1 << ")^T\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "TT," << tk+1 << "," << ti << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif

						SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk+1,tk+1), A(ti,tk+1) );
						#ifdef COUT
						#pragma omp critical
						cout << "SSRFB_L( (" << tk+1 << "," << tk+1 << "," << tk << "), (" << ti << "," << tk+1 << "," << tk << ") )\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "SL," << ti << "," << tk+1 << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					// (f), (g)
					for (int tj=tk+2; tj<ti; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) depend(inout: Ap[tj][tk+1], Ap[ti][tj])
						{
							BMatrix tmp0(NB,NB,IB);

							tileTransCopy(A(tj,tk+1),tmp0);
							SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), &tmp0, A(ti,tj) );
							tileTransCopy(tmp0,A(tj,tk+1));
							#ifdef COUT
							#pragma omp critical
							cout << "SSRFB_L( (" << tj << "," << tk+1 << "," << tk << ")^T, (" << ti << "," << tj << "," << tk << ") )\n";
							#endif
							#ifdef ANIM
							#pragma omp critical
							cout << "SL," << ti << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
							#endif
						}
					}

					//////////////////////////////////
					// (h)
					#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) depend(inout: tmp, Ap[ti][ti])
					{
						SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), &tmp, A(ti,ti) );
						#ifdef COUT
						#pragma omp critical
						cout << "SSRFB_L( T(" << tk+1 << "," << ti << "," << tk << "), (" << ti << "," << ti << "," << tk << ") )\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "SL," << ti << "," << ti << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					//////////////////////////////////
					// (i)
					#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) depend(inout: Ap[tk+1][tk+1], tmp)
					{
						SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tk+1,tk+1), &tmp );
						#ifdef COUT
						#pragma omp critical
						cout << "SSRFB_R( (" << tk+1 << "," << tk+1 << "," << tk << "), T(" << tk+1 << "," << ti << "," << tk << ") )\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "SR," << tk+1 << "," << ti << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					// (k)
					for (int tj=ti; tj < NT; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) depend(inout: Ap[tj][tk+1], Ap[tj][ti])
						{
							SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tj,tk+1), A(tj,ti) );
							#ifdef COUT
							#pragma omp critical
							cout << "SSRFB_R( (" << tj << "," << tk+1 << "," << tk << "), (" << tj << "," << ti << "," << tk << ") )\n";
							#endif
							#ifdef ANIM
							#pragma omp critical
							cout << "SR," << tj << "," << ti << "," << tk << "," << omp_get_wtime() - ttime << endl;
							#endif
						}
					} // j-LOOP END
				} // i-LOOP END
			} // k-LOOP END
		} // single section END
	} // parallel section END
}
