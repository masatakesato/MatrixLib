#ifndef HOUSEHOLDER_H
#define	HOUSEHOLDER_H


#include	"MatrixLib.h"




template< typename T >
void Householder( IMatrix<T>& a, const IMatrix<T>& v )
{
	T v2 = 0;
	// calc vtv
	for( int i=0; i<v.numRows(); ++i)
		v2 += pow( v(i, 0), 2 );

	T coeff = 2 / v2;
	for( int i=0; i<v.numRows(); ++i )
	{
		for( int j=i; j<v.numRows(); ++j )
		{
			// I - 2 * v(i, j) * v(j, i) / vtv
			a(i, j) = T(i==j) - coeff * v(i, 0) * v(j, 0);
			a(j, i) = a(i, j);
		}

	}


}

// http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?%B8%C7%CD%AD%C3%CD/%B8%C7%CD%AD%A5%D9%A5%AF%A5%C8%A5%EB

#endif	// !HOUSEHOLDER_H //