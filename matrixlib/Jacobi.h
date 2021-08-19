#ifndef JACOBI_H
#define	JACOBI_H


//#include	"IMatrix.h"
//#include	"DynamicMatrix.h"
#include	"MatrixLib.h"


template< typename T >
void Jacobi( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	/*IMatrix<T>*/DynamicMatrix<T> x_[] = { x, x };
	unsigned int curr = 0;
	DynamicMatrix<T> &x_curr = x_[curr], &x_next = x_[curr^1];


	for( int k=0; k<maxiter; k++ )
	{
		for( int i=0; i<x.numRows(); i++ )
		{
			// Calculate x_k+1
			T& x_i = x_[ curr^1 ](i, 0);
			x_i = b( i, 0 );
			for( int j=0; j<A.numCols(); j++ )
			{
				if( i==j ) continue;
				x_i -= A(i, j) * x_curr(j, 0);
			}
			x_i /= A(i, i);
			
		}// end of i loop
		
		// Check convergence
		T err = 0;
		for( int i=0; i<x.numRows(); i++ )
			err += abs( x_next(i, 0) - x_curr(i, 0) );

		if( err < eps )	break;

		curr ^= 1;
		x_next = x_[ curr^1 ];
		x_curr = x_[ curr ];

	}// end of k loop

	x.Copy(x_next);

}



#endif // !JACOBI_H
