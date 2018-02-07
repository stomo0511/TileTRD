//  RightLookingTRD_Task
//
//  Created by T. Suzuki on 2018/01/30.
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

using namespace std;

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
				//////////////////////////////////////////////////////////////
				// Left
				#pragma omp task depend(inout:Ap[tk+1][tk]) \
								 depend(out:Tp[tk+1][tk])
				{
					GEQRT( A(tk+1,tk), T(tk+1,tk) );
					#ifdef COUT
					#pragma omp critical
					cout << "GEQRT_L(" << tk+1 << "," << tk << "," << tk << ")\n";
					#endif
					#ifdef ANIM
					cout << "GL," << tk+1 << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

				for  (int tj=tk+1; tj < NT; tj++)
				{
					#pragma omp task depend(in:Ap[tk+1][tk], Tp[tk+1][tk]) \
					                 depend(inout:Ap[tk+1][tj])
					{
						LARFB( PlasmaLeft, PlasmaTrans,   A(tk+1,tk), T(tk+1,tk), A(tk+1,tj) );
						#ifdef COUT
						#pragma omp critical
						cout << "LARFB_L(" << tk+1 << "," << tj << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "LL," << tk+1 << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}
				} // j-LOOP END

				for (int ti=tk+2; ti < MT; ti++)
				{
					#pragma omp task depend(inout:Ap[tk+1][tk], Ap[ti][tk]) \
					                 depend(out:Tp[ti][tk])
					{
						TSQRT( A(tk+1,tk), A(ti,tk), T(ti,tk) );
						#ifdef COUT
						#pragma omp critical
						cout << "TSQRT_L( A(" << tk+1 << "," << tk << "), A(" << ti << "," << tk << ") )\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "TL," << ti << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					for (int tj=tk+1; tj < NT; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) \
						                 depend(inout: Ap[tk+1][tj], Ap[ti][tj])
						{
							SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk+1,tj), A(ti,tj) );
							#ifdef COUT
							#pragma omp critical
							cout << "SSRFB_L( A(" << tk+1 << "," << tj << "), A(" << ti << "," << tj << ") )\n";
							#endif
							#ifdef ANIM
							#pragma omp critical
							cout << "SL," << ti << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
							#endif
						}
					} // j-LOOP END
				} // i-LOOP END

				//////////////////////////////////////////////////////////////
				// Right
				#pragma omp task depend(in:Ap[tk+1][tk])
				{
					// Copy A(tk+1,tk)^T to A(tk,tk+1)
					for (int i=0; i<A(tk+1,tk)->m(); i++)
						for (int j=0; j<A(tk+1,tk)->n(); j++) {
							A(tk,tk+1)->Set_Val( j, i, A(tk+1,tk)->Get_Val(i,j) );
						}
					#ifdef COUT
					#pragma omp critical
					cout << "GEQRT_R(" << tk << "," << tk+1 << "," << tk << ")\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "GR," << tk << "," << tk+1 << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

				for (int tj=tk+1; tj < NT; tj++)
				{
					#pragma omp task depend(in:Ap[tk+1][tk], Tp[tk+1][tk]) \
					                 depend(inout:Ap[tj][tk+1])
					{
						LARFB( PlasmaRight, PlasmaNoTrans, A(tk+1,tk), T(tk+1,tk), A(tj,tk+1) );
						#ifdef COUT
						#pragma omp critical
						cout << "LARFB_R(" << tj << "," << tk+1 << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "LR," << tj << "," << tk+1 << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}
				} // j-LOOP END

				for (int ti=tk+2; ti < MT; ti++)
				{
					#pragma omp task depend(in:Ap[tk+1][tk], Ap[ti][tk])
					{
						// Copy A(tk+1,tk)^T to A(tk,tk+1)
						for (int i=0; i<A(tk+1,tk)->m(); i++)
							for (int j=0; j<A(tk+1,tk)->n(); j++) {
								A(tk,tk+1)->Set_Val( j, i, A(tk+1,tk)->Get_Val(i,j) );
							}

						// Copy A(ti,tk)^T to A(tk,ti)
						for (int i=0; i<A(ti,tk)->m(); i++)
							for (int j=0; j<A(ti,tk)->n(); j++) {
								A(tk,ti)->Set_Val( j, i, A(ti,tk)->Get_Val(i,j) );
							}
						#ifdef COUT
						#pragma omp critical
						cout << "TSQRT_R(" << tk << "," << ti << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "TR," << tk << "," << ti << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					for (int tj=tk+1; tj < NT; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk], Tp[ti][tk]) \
						                 depend(inout: Ap[tj][tk+1], Ap[tj][ti])
						{
							SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tj,tk+1), A(tj,ti) );
							#ifdef COUT
							#pragma omp critical
							cout << "SSRFB_R( A(" << tj << "," << tk+1 << "), A(" << tj << "," << ti << ") )\n";
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
