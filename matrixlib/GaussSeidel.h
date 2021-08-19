#ifndef GAUSS_SEIDEL_H
#define	GAUSS_SEIDEL_H


#include	"MatrixLib.h"




template< typename T >
void GaussSeidel( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{

	for( int k=0; k<maxiter; k++ )
	{
		T err = 0;

		for( int i=0; i<x.numRows(); i++ )
		{
			// Calculate x_k+1
			T& x_i = x(i, 0);
			T tmp = x_i;

			x_i = b( i, 0 );
			for( int j=0; j<A.numCols(); j++ )
			{
				if( i==j ) continue;
				x_i -= A(i, j) * x(j, 0);
			}
			x_i /= A(i, i);

			err += abs( tmp - x_i );

		}// end of i loop

		 // Check convergence
		if( err < eps )	break;

	}// end of k loop

}






#endif // !GAUSS_SEIDEL_H
