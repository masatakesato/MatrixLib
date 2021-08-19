#ifndef GRAM_SCHMIDT_H
#define	GRAM_SCHMIDT_H


#include	"MatrixLib.h"


// http://zakii.la.coocan.jp/signal/35_qr_gram_schmidt.htm
// https://www.math.uci.edu/~ttrogdon/105A/html/Lecture23.html



// Modified Gram-Schmidt orthogonalization
template< typename T >
void ModifiedGramSchmidt( IMatrix<T>& v, IMatrix<T>&q, const IMatrix<T>& x )
{
	int n = x.numCols();// num of vectors
	int m = x.numRows();// vector dimension

	for( int j=0; j<n; ++j )
		for( int l=0; l<m; ++l ) v(l, j) = x(l, j);

	for( int j=0; j<n; ++j )
	{
		// qj = vj / ||vj||2
		T vj_len = 0;
		for( int l=0; l<m; ++l )
			vj_len += pow( v(l, j), 2 );
		vj_len = sqrt( vj_len );
		for( int l=0; l<m; ++l )
			q(l, j) = v(l, j) / vj_len;

		for( int k=j+1; k<n; ++k )
		{
			// vk = vk - (qj,  vk) * qj
			T dot_qj_vk = 0;
			for( int l=0; l<m; ++l )
				dot_qj_vk += q(l, j) * v(l, k);
			for( int l=0; l<m; ++l )
				v(l, k) -= dot_qj_vk * q(l, j);
		}// end of k loop

	}// end of i loop

}



#endif // !GRAM_SCHMIDT_H
