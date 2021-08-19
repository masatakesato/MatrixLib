#ifndef SUBSTITUTION_H
#define	SUBSTITUTION_H



#include	"IMatrix.h"



// Solve Lx = b
//   L: lower triangular matrix
//   x: unknown variables
//   b: solutions
template< typename T >
inline void ForwardSubstitution( const IMatrix<T>& L, IMatrix<T>& x, const IMatrix<T>& b )
{
	for( int k=0; k<b.numCols(); ++k )
	{
		for( int i=0; i<L.numRows(); ++i )
		{
			x(i, k) = b(i, k);
			for( int j=0; j<i; ++j )
				x(i, k) -= L(i, j) * x(j, k);
			x(i, k) /= L(i, i);
		}
	}
}


// solves column k only
//template< typename T >
//inline void ForwardSubstitution( const Matrix<T>& L, Matrix<T>& x, const Matrix<T>& b, int k )
//{
//	for( int i=0; i<L.numRows(); ++i )
//	{
//		x(i, k) = b(i, k);
//		for( int j=0; j<i; ++j )
//			x(i, k) -= L(i, j) * x(j, k);
//		x(i, k) /= L(i, i);
//	}
//}



// Solve Ux = b
//   U: upper triangular matrix
//   x: unknown variables
//   b: solutions
template< typename T >
inline void BackwardSubstitution( const IMatrix<T>& U, IMatrix<T>& x, const IMatrix<T>& b )
{
	for( int k=0; k<b.numCols(); ++k )
	{
		for( int i=U.numRows()-1; i>=0; --i )
		{
			x(i, k) = b(i, k);
			for( int j=i+1; j<U.numCols(); ++j )
				x(i, k) -= U(i, j) * x(j, k);
			x(i, k) /= U(i, i);
		}
	}
}


// solves column k only
//template< typename T >
//inline void BackwardSubstitution( const Matrix<T>& U, Matrix<T>& x, const Matrix<T>& b, int k )
//{
//	for( int i=U.numRows()-1; i>=0; --i )
//	{
//		x(i, k) = b(i, k);
//		for( int j=i+1; j<U.numCols(); ++j )
//			x(i, k) -= U(i, j) * x(j, k);
//		x(i, k) /= U(i, i);
//	}
//}



#endif // !SUBSTITUTION_H