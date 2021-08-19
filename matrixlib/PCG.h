#ifndef PRECONDITIONED_CONJUGATE_GRADIENT_H
#define	PRECONDITIONED_CONJUGATE_GRADIENT_H


//#include	"Matrix.h"
#include	"Solver.h"




template< typename T >
void PCG( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, Solver<T>& precond, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				z( x.numRows(), x.numCols() ),
				z1( x.numRows(), x.numCols() ),
				p( x.numRows(), x.numCols() ),
				Ap( x.numRows(), x.numCols() );


	//============== Initialize r and p =============//

	// r0 = b - Ax0
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	
	}

	// z0 = (LDLt)^-1 * r
	precond.Solve( z, r );

	// p0=z0
	p = z;

	T rz = Dot(r, z), rz1;


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// A * p
		Multiply( Ap, A, p );
		

		// alpha
		T alpha = rz / Dot(p, Ap);
		
		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
		}

		// convergence check
		if( sqrt( Dot(r, r) ) < eps )
			break;

		// z_k+1
		precond.Solve( z1, r );

		// beta
		rz1 = Dot(r, z1);
		T beta = rz1 / rz;

		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = z1(i, 0) + beta * p(i, 0);

		rz = rz1;

	}// end of k loop

}



// Preconditioned BiCG. http://ir.lib.u-ryukyu.ac.jp/bitstream/20.500.12000/2218/1/No55p17.pdf
template< typename T >
void PBiCG( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, Solver<T>& precond, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				z1( x.numRows(), x.numCols() ),
				z2( x.numRows(), x.numCols() ),
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
	}

	precond.Solve( p, r );//p0 = K^-1 * r0
	ps = p;// ps0 = (LDLt)^-1 * rs0

	T	rsz = Dot(rs, p),
		rsz1;

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		Multiply( Ap, A, p );
		Multiply( Atps, At, ps );
		
		T alpha = rsz / Dot(Ap, ps);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
			rs(i, 0) -= alpha * Atps(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )
			break;

		precond.Solve( z1, r );// z1 = K^-1 * r_k+1
		rsz1 = Dot(rs, z1);
		T beta = rsz1 / rsz;

		precond.Solve( z2, rs );// z2 = r*_k+1 * K^-1

		for( int i=0; i<p.numRows(); i++ )
		{
			p(i, 0) = z1(i, 0) + beta * p(i, 0);
			ps(i, 0) = z2(i, 0) + beta * ps(i, 0);
		}

		rsz = rsz1;
	}// end of k loop

}




template< typename T >
void PBiCGSTAB( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, Solver<T>& precond, int maxiter=100, double eps=1.0e-8 )
{
	DynamicMatrix<T>	r( x.numRows(), x.numCols() ),
				Ks( x.numRows(), x.numCols() ),
				p( x.numRows(), x.numCols() ),
				Kp( x.numRows(), x.numCols() ),
				AKp( x.numRows(), x.numCols() ),
				AKs( x.numRows(), x.numCols() );

	DynamicMatrix<T>	rs( x.numRows(), x.numCols() ),
				s( x.numRows(), x.numCols() );


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		rs(i, 0) = r(i, 0);		// r* = r

		p(i, 0) = /*1;*/r(i, 0);		// p = r
	}

	T	rrs = Dot(r, rs),
		rrs1;

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// Kp = M * p -> M^-1 * Kp = p
		precond.Solve( Kp, p );
		Multiply( AKp, A, Kp );

		T alpha = rrs / Dot(AKp, rs);

		for( int i=0; i<s.numRows(); i++ )
			s(i, 0) = r(i, 0) - alpha * AKp(i, 0);

		precond.Solve( Ks, s );
		Multiply( AKs, A, Ks );

		T aks_len = Dot(AKs, AKs);
		T w = aks_len==0 ? 1 : Dot(AKs, s) / aks_len;
		//T w = Dot(AKs, s) / Dot(AKs, AKs);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * Kp(i, 0) + w * Ks(i, 0);
			r(i, 0) = s(i, 0) - w * AKs(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )
			break;

		rrs1 = Dot(r, rs);
		T beta = rrs1 / rrs * alpha / w;
		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = r(i, 0) + beta * ( p(i, 0) - w * AKp(i, 0) );

		rrs = rrs1;

	}// end of k loop

}









//##############################################################################################################//
//							Deprecated: Hard-coded preconditioner version. 2021.03.04							//
//##############################################################################################################//

/*

#include	"Cholesky.h"


template< typename T >
void SolveTridiagonal( const Matrix<T>& Lower, const Matrix<T>& Upper, Matrix<T>& x, const Matrix<T>& b )
{
// Solve Lower * y = b 
Matrix<T> y( x.numRows(), x.numCols() );
ForwardSubstitution( Lower, y, b );

// Solve Upper * x = y
BackwardSubstitution( Upper, x, y );
}



// https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

// Incomplete cholesky decomposition + CG
template< typename T >
void ICCG( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	Matrix<T>	r( x.numRows(), x.numCols() ),
		z( x.numRows(), x.numCols() ),
		z1( x.numRows(), x.numCols() ),
		p( x.numRows(), x.numCols() ),
		Ap( x.numRows(), x.numCols() );

	Matrix<T>	L( A.numRows(), A.numCols() ),
		D( A.numRows(), A.numCols() ),
		LD( A.numRows(), A.numCols() ),
		Lt( A.numRows(), A.numCols() );


	//============== Initialize r and p =============//

	// r0 = b - Ax0
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	
	}

	// z0 = (LDLt)^-1 * r
	IncompleteCholesky( L, D, A );
	Multiply( LD, L, D );
	Transpose( Lt, L );

	SolveTridiagonal( LD, Lt, z, r );
	// Solve L*y=r (y=D*Lt*p) and get variavle y.
	// Solve D*Lt*p=y and get p

	// p0=z0
	p = z;

	T rz = Dot(r, z), rz1;


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// A * p
		Multiply( Ap, A, p );

		// alpha
		T alpha = rz / Dot(p, Ap);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
		}

		// convergence check
		if( sqrt( Dot(r, r) ) < eps )
			break;

		// z_k+1
		SolveTridiagonal( LD, Lt, z1, r );

		// beta
		rz1 = Dot(r, z1);
		T beta = rz1 / rz;

		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = z1(i, 0) + beta * p(i, 0);

		rz = rz1;

	}// end of k loop

}



// https://www.researchgate.net/figure/The-Preconditioned-Bi-Conjugate-Gradient-Stabilized-Method_fig17_267701113
// https://www.jstage.jst.go.jp/article/jsiamt/23/2/23_KJ00008726346/_pdf
// https://en.wikipedia.org/wiki/Biconjugate_gradient_method
// https://aip.scitation.org/doi/pdf/10.1063/1.4823087

template< typename T >
void ICBiCG( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	Matrix<T>	r( x.numRows(), x.numCols() ),
		z1( x.numRows(), x.numCols() ),
		z2( x.numRows(), x.numCols() ),
		p( x.numRows(), x.numCols() ),
		Ap( x.numRows(), x.numCols() );

	Matrix<T>	rs( x.numRows(), x.numCols() ),
		ps( x.numRows(), x.numCols() ),
		Atps( x.numRows(), x.numCols() );

	Matrix<T>	L( A.numRows(), A.numCols() ),
		D( A.numRows(), A.numCols() ),
		Lt( A.numRows(), A.numCols() ),
		LD( A.numRows(), A.numCols() );

	Matrix<T>	At = A;
	At.Transpose();

	IncompleteCholesky(L, D, A);
	Multiply( LD, L, D );
	Transpose( Lt, L );


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		rs(i, 0) = r(i, 0);		// r* = r
	}

	SolveTridiagonal( LD, Lt, p, r );//p0 = K^-1 * r0
	ps = p;// ps0 = (LDLt)^-1 * rs0

	T	rsz = Dot(rs, p),
		rsz1;

	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		Multiply( Ap, A, p );
		Multiply( Atps, At, ps );

		T alpha = rsz / Dot(Ap, ps);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * p(i, 0);
			r(i, 0) -= alpha * Ap(i, 0);
			rs(i, 0) -= alpha * Atps(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )
			break;

		SolveTridiagonal( LD, Lt, z1, r );// z1 = K^-1 * r_k+1
		rsz1 = Dot(rs, z1);
		T beta = rsz1 / rsz;

		SolveTridiagonal( LD, Lt, z2, rs );// z2 = r*_k+1 * K^-1

		for( int i=0; i<p.numRows(); i++ )
		{
			p(i, 0) = z1(i, 0) + beta * p(i, 0);
			ps(i, 0) = z2(i, 0) + beta * ps(i, 0);
		}

		rsz = rsz1;
	}// end of k loop

}




// http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95

template< typename T >
void ICBiCGSTAB( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b, int maxiter=100, double eps=1.0e-8 )
{
	Matrix<T>	r( x.numRows(), x.numCols() ),
		Ks( x.numRows(), x.numCols() ),
		p( x.numRows(), x.numCols() ),
		Kp( x.numRows(), x.numCols() ),
		AKp( x.numRows(), x.numCols() ),
		AKs( x.numRows(), x.numCols() );

	Matrix<T>	rs( x.numRows(), x.numCols() ),
		s( x.numRows(), x.numCols() );

	Matrix<T>	L( A.numRows(), A.numCols() ),
		D( A.numRows(), A.numCols() ),
		LD( A.numRows(), A.numCols() ),
		Lt( A.numRows(), A.numCols() );


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); i++ )
	{
		T ax = 0;
		for( int j=0; j<A.numCols(); j++ )
			ax += A(i, j) * x(j, 0);

		r(i, 0) = b(i, 0) - ax;	// r = b - Ax
		rs(i, 0) = r(i, 0);		// r* = r

		p(i, 0) = 1;//r(i, 0);		// p = r
	}

	T	rrs = Dot(r, rs),
		rrs1;

	IncompleteCholesky(L, D, A);
	Multiply( LD, L, D );
	Transpose( Lt, L );


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		SolveTridiagonal( LD, Lt, Kp, p );
		Multiply( AKp, A, Kp );

		T alpha = rrs / Dot(AKp, rs);

		for( int i=0; i<s.numRows(); i++ )
			s(i, 0) = r(i, 0) - alpha * AKp(i, 0);

		SolveTridiagonal( LD, Lt, Ks, s );
		Multiply( AKs, A, Ks );

		T aks_len = Dot(AKs, AKs);
		T w = aks_len==0 ? 1 : Dot(AKs, s) / aks_len;// T w = Dot(AKs, s) / Dot(AKs, AKs);

		for( int i=0; i<x.numRows(); i++ )
		{
			x(i, 0) += alpha * Kp(i, 0) + w * Ks(i, 0);
			r(i, 0) = s(i, 0) - w * AKs(i, 0);
		}

		if( sqrt( Dot(r, r) ) < eps )
			break;

		rrs1 = Dot(r, rs);
		T beta = rrs1 / rrs * alpha / w;
		for( int i=0; i<p.numRows(); i++ )
			p(i, 0) = r(i, 0) + beta * ( p(i, 0) - w * AKp(i, 0) );

		rrs = rrs1;

	}// end of k loop

}

*/


#endif // !PRECONDITIONED_CONJUGATE_GRADIENT_H
