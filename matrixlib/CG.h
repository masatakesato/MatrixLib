#ifndef CONJUGATE_GRADIENT_H
#define	CONJUGATE_GRADIENT_H


#include	"MatrixLib.h"




template< typename T >
void CG( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				p( x.numRows(), x.numCols() ),
				Ap( x.numRows(), x.numCols() );

	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		p(i, 0) = r(i, 0);		// p = r
	}

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// A * p
		Multiply( Ap, A, p );

		T rr = Dot(r, r);
		T alpha = rr / Dot(p, Ap);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
		}

		T rr_next = Dot(r, r);

		if( sqrt( rr_next ) < eps )	break;

		T beta = rr_next / rr;

		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = r(i, 0) + beta * p(i, 0);

	}// end of k loop

}



//template< typename T >
//void CG( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b, int maxiter=100, double eps=1.0e-8 )
//{
//	Matrix<T> r = b - A * x, p = r;
//
//	for( int k=0; k<maxiter; k++ )
//	{
//		auto y = A * p;
//
//		T rr = Dot(r, r);
//		T alpha = rr / Dot(p, y);
//
//		x += alpha * p;
//		r -= alpha * y;
//
//		T rr_next = Dot(r, r);
//		
//		if( sqrt( rr_next ) < eps )	break;
//
//		auto beta = rr_next / rr;
//
//		p = r + beta * p;
//
//	}// end of k loop
//
//}





template< typename T >
void BiCG( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				p( x.numRows(), x.numCols() ),
				Ap( x.numRows(), x.numCols() );

	DynamicMatrix<T>	rs( x.numRows(), x.numCols() ),
				ps( x.numRows(), x.numCols() ),
				Atps( x.numRows(), x.numCols() );

	DynamicMatrix<T>	At = A;
	At.Transpose();


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		rs(i, 0) = r(i, 0);		// r* = r

		p(i, 0) = r(i, 0);		// p = r
		ps(i, 0) = rs(i, 0);	// p* = r
	}

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		Multiply( Ap, A, p );
		Multiply( Atps, At, ps );

		T rrs = Dot(r, rs);
		T alpha = rrs / Dot(Ap, ps);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
			rs(i, 0) -= alpha * Atps(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )	break;

		T beta = Dot(r, rs) / rrs;

		for( int i=0; i<p.numRows(); i++ )
		{
			p(i, 0) = r(i, 0) + beta * p(i, 0);
			ps(i, 0) = rs(i, 0) + beta * ps(i, 0);
		}

	}// end of k loop

}




// http://www-section.cocolog-nifty.com/blog/2008/12/bicgstab-30b8-1.html
// http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
// https://qiita.com/fukuroder/items/4b708524783192fc2018

template< typename T >
void BiCGSTAB( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				p( x.numRows(), x.numCols() ),
				Ap( x.numRows(), x.numCols() );

	DynamicMatrix<T>	rs( x.numRows(), x.numCols() ),
				s( x.numRows(), x.numCols() ),
				As( x.numRows(), x.numCols() );


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		rs(i, 0) = r(i, 0);		// r* = r

		p(i, 0) = r(i, 0);		// p = r
	}

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		Multiply( Ap, A, p );

		T rrs = Dot(r, rs);
		T alpha = rrs / Dot(Ap, rs);
		
		for( int i=0; i<s.numRows(); i++ )
			s(i, 0) = r(i, 0) - alpha * Ap(i, 0);

		Multiply( As, A, s );

		T w = Dot(As, s) / Dot(As, As);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0) + w * s(i, 0);
			r(i, 0) = s(i, 0) - w * As(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )	break;

		T beta = Dot(r, rs) / rrs * alpha / w;
		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = r(i, 0) + beta * ( p(i, 0) - w * Ap(i, 0) );

	}// end of k loop

}




#endif // !CONJUGATE_GRADIENT_H
