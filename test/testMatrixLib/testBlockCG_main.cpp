#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/CG.h>
#include	<matrixlib/PCG.h>

#include	<matrixlib/LU.h>
#include	<matrixlib/QR.h>
#include	<matrixlib/Cholesky.h>
#include	<matrixlib/SOR.h>



template< typename T >
inline T FrobeniusNorm( const IMatrix<T>& in )
{
	T sumSqrd = 0;
	auto pin = in.pDataAs<T>();
	for( int i=0; i<in.numElms(); ++i )
		sumSqrd += pow( pin[i], 2 );

	return sqrt(sumSqrd);
}




// https://row1.ca/pdf/dc-resistivity-block-cg.pdf
// https://github.com/lkeegan/blockCG/blob/master/inc/block_solvers.hpp


// ランク落ち対策
// Initial Deflation: Bの代わりにQR分解したQで反復計算. 最後にRをかける
// A * X = B, B = Q * R
// X = A^-1 * Q * R
// y = A^-1 * Q
// X = y * R

// Internal Deflation: 反復計算中はPとRのランクを監視する.



// 最初からBを直交基底に置き換えて計算する


template< typename T >
void BlockCG( const IMatrix<T>& A, IMatrix<T>& X, const IMatrix<T>& B, int maxiter=100, double eps=1.0e-8 )
{
	// M: row/col sizeof A
	// N: col size of B
	DynamicMatrix<T>	R( X.numRows(), X.numCols() ),// M x N
				P( X.numRows(), X.numCols() ),// M x N
				AP( X.numRows(), X.numCols() ),// M x N

				Alpha( X.numCols(), X.numCols() ),// N x N
				Beta( X.numCols(), X.numCols() );// N x N

	DynamicMatrix<T>	RtR( X.numCols(), X.numCols() ),// N x N
				RtR_next( X.numCols(), X.numCols() );// N x N

	DynamicMatrix<T>	bufNxN( X.numCols(), X.numCols() ),
				bufMxN( X.numRows(), X.numCols() );

	PivotLU_Solver<T>	solver(Alpha, 2);
	

	//============== Initialize R and P =============//
	for( int i=0; i<R.numRows(); ++i )
	{
		for( int j=0; j<R.numCols(); ++j )
		{
			T ax = 0;
			for( int k=0; k<A.numCols(); ++k )
				ax += A(i, k) * X(k, j);

			R(i, j) = B(i, j) - ax;	// R0 = B - AX
			P(i, j) = R(i, j);		// P0 = R0
		}
	}
	
	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// AP = A * P
		Multiply( AP, A, P );

		// Solve (P^t * AP) * Alpha = R^t * R
		//T rr = Dot(r, r);
		//T alpha = rr / Dot(p, Ap);
		TransposeMultiply( bufNxN, P, AP );
		solver.SetCoeff( bufNxN );
		TransposeMultiply( RtR, R );// R^t * R
		solver.Solve( Alpha, RtR );

		// X_next = X_current + P*Alpha
		Multiply( bufMxN, P, Alpha );
		Add( X, bufMxN );

		// R_next = R_current - AP*Alpha
		Multiply( bufMxN, AP, Alpha );
		Subtract( R, bufMxN );

		// Check conversion
		//T rr_next = Dot(r, r);
		//if( sqrt( rr_next ) < eps )	break;
		TransposeMultiply( RtR_next, R );// R_next^t * R_next
		if( FrobeniusNorm( RtR_next ) < eps ) break;

		// Solve (R^t * R) * Beta = R_next^t * R_next
		//T beta = rr_next / rr;
		solver.SetCoeff( RtR );
		solver.Solve( Beta, RtR_next );

		// P_next = R_next + P * Beta;
		//for( int i=0; i<p.numRows(); i++ )
		//	p(i, 0) = r(i, 0) + beta * p(i, 0);
		Multiply( bufMxN, P, Beta );
		Add( P, R, bufMxN );

	}// end of k loop

}



// https://www.cs.odu.edu/~yaohang/portfolio/BIT2017.pdf

template< typename T >
void BreakdownFreeBlockCG( const IMatrix<T>& A, IMatrix<T>& X, const IMatrix<T>& B, Solver<T>& precond, int maxiter=100, double eps=1.0e-9 )
{
	DynamicMatrix<T>	R( X.numRows(), X.numCols() ),// M x N
				P( X.numRows(), X.numCols() ),// M x N
				AP( X.numRows(), X.numCols() ),// M x N

				Alpha( X.numCols(), X.numCols() ),// N x N
				BetaNeg( X.numCols(), X.numCols() ),// N x N. -Beta

				Z( X.numRows(), X.numCols() );// M x N
				
	DynamicMatrix<T>	bufNxN( X.numCols(), X.numCols() ),
				bufMxN( X.numRows(), X.numCols() );


	QRColumnPivotHouseholder_Solver<T> orth(Z);
	PivotLU_Solver<T>	solver(Alpha, 2);
				

	//============== Initialize R0, Z0 and P0 =============//
	// R0 = B - AX0
	for( int i=0; i<R.numRows(); ++i )
	{
		for( int j=0; j<R.numCols(); ++j )
		{
			T ax = 0;
			for( int k=0; k<A.numCols(); ++k )
				ax += A(i, k) * X(k, j);

			R(i, j) = B(i, j) - ax;	// R0 = B - AX
		}
	}

	// Z0 = M * R
	precond.Solve(Z, R);

	// P = orth(Z0)
	orth.Init(Z);
	for( int i=0; i<P.numCols(); ++i )
		P.CopyColumn( orth.Q(), i, i );

	//================== Iterate ===================//
	for( int k=0; k<maxiter; ++k )
	{
		// AP = A * P
		Multiply( AP, A, P );

		// Solve (P^t * AP) * Alpha = P^t * R
		TransposeMultiply( bufNxN, P, AP );
		solver.SetCoeff( bufNxN );
		TransposeMultiply( bufNxN, P, R );// P^t * R
		solver.Solve( Alpha, bufNxN );
		
		// X_next = X_current + P*Alpha
		Multiply( bufMxN, P, Alpha );
		Add( X, bufMxN );
		
		// R_next = R_current - AP*Alpha
		Multiply( bufMxN, AP, Alpha );
		Subtract( R, bufMxN );

		// Check conversion
		TransposeMultiply( bufNxN, R );// R_next^t * R_next
		if( FrobeniusNorm( bufNxN ) < eps )
			break;
		
		// Z_next = M * R_next
		precond.Solve( Z, R );
		
		// Solve (P^t * AP) * Beta = -Q^t * Z_next
		TransposeMultiply( bufNxN, AP, Z );// -Q^t * Z_next
		solver.Solve( BetaNeg, bufNxN );// get NEGATIVE Beta!
		
		// P = orth(Z_next + P*Beta) -> P = orth(Z_next - P * BetaNeg )
		Multiply( bufMxN, P, BetaNeg );// P * NEGATIVE Beta
		Subtract( Z, bufMxN );// Z - (P * NEGATIVE Beta)

		orth.Init(Z);
		for( int i=0; i<P.numCols(); ++i )
			P.CopyColumn( orth.Q(), i, i );

	}// end of k loop

}




// https://github.com/nmoteki/block-Krylov-linear-solvers

template< typename T >
void BlockBiCGSTAB( const IMatrix<T>& A, IMatrix<T>& X, const IMatrix<T>& B, int maxiter=100, double eps=1.0e-9 )
{

	DynamicMatrix<T>	R( X.numRows(), X.numCols() ),// M x N
				P( X.numRows(), X.numCols() ),// M x N
				AP( X.numRows(), X.numCols() ),// M x N

				Alpha( X.numCols(), X.numCols() ),// N x N
				BetaNeg( X.numCols(), X.numCols() );// N x N. -Beta

	DynamicMatrix<T>	Rs( X.numRows(), X.numCols() ),// M x N
				S( X.numRows(), X.numCols() ),// M x N
				AS( X.numRows(), X.numCols() ),// M x N
				RsAP( X.numCols(), X.numCols() );// N x N. Rs^t * AP

	DynamicMatrix<T>	bufNxN( X.numCols(), X.numCols() ),// N x N
				bufMxN( X.numRows(), X.numCols() );// M x N

	PivotLU_Solver<T>	solver(Alpha, 2);


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); ++i )
	{
		for( int j=0; j<R.numCols(); ++j )
		{
			T ax = 0;
			for( int k=0; k<A.numCols(); ++k )
				ax += A(i, k) * X(k, j);

			R(i, j) = B(i, j) - ax;	// R0 = B - AX0
			Rs(i, j) = R(i, j);	// R* = R

			P(i, j) = R(i, j);		// P0 = R0
		}
	}


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// AP = A * Pk
		Multiply( AP, A, P );

		// Solve (Rs^t * AP) * Alpha = Rs^t * R
		TransposeMultiply( RsAP, Rs, AP );
		solver.SetCoeff( RsAP );
		TransposeMultiply( bufNxN, Rs, R );// Rs^t * R
		solver.Solve( Alpha, bufNxN );

		// Sk = Rk - AP * Alpha
		Multiply( bufMxN, AP, Alpha );
		Subtract( S, R, bufMxN );

		// AS = A * Sk
		Multiply( AS, A, S );

		// z = Tr( AS^t * S ) / Tr( AS^t * AS )
		TransposeMultiply( bufNxN, AS, S );
		T z = Tr( bufNxN );// z = Tr( AS^t * S )
		TransposeMultiply( bufNxN, AS );
		T tt = Tr( bufNxN );
		z /= tt==0 ? 1.0e-9 : tt;// z /= Tr( AS^t * AS )

		// X_next += Pk * Alpha + z * Sk
		Multiply( bufMxN, P, Alpha );// Pk * Alpha
		Add( X, bufMxN );// X_next += Pk * Alpha
		Scale( bufMxN, z, S );//z * Sk
		Add( X, bufMxN );// X_next += z * Sk

		// R_next = Sk - z * ASk
		Scale( bufMxN, z, AS );// z * ASk
		Subtract( R, S, bufMxN );//Sk - z * ASk

		// Check conversion
		TransposeMultiply( bufNxN, R );// R_next^t * R_next
		if( FrobeniusNorm( bufNxN ) < eps )
			break;

		// Solve (Rs^t * AP) * Beta = -Rs^t * ASk
		TransposeMultiply( bufNxN, Rs, AS );// Rs^t * ASk
		solver.Solve( BetaNeg, bufNxN );// Get NEGATIVE Beta!

		// P_next = R_next + (P_current - z * APk) * Beta   ->   P_next = R_next + (-P_current + z * APk) * BetaNeg
		AddScaled( bufMxN, P, (T)-1, AP, z );// (-P_current + z * APk)
		Multiply( P, bufMxN, BetaNeg );// P_next = (-P_current + z * APk) * BetaNeg
		Add( P, R );// P_next += R_next

	}// end of k loop

}




// http://computics-material.jp/jpn/symposium/20110301/pdf/05_1450_1500_tadano.pdf // Block version with QR decomposition

template< typename T >
void BlockBiCGSTAB_QR( const IMatrix<T>& A, IMatrix<T>& X, const IMatrix<T>& B, int maxiter=100, double eps=1.0e-9 )
{
	DynamicMatrix<T>	R( X.numRows(), X.numCols() ),// M x N
				P( X.numRows(), X.numCols() ),// M x N
				AQ( X.numRows(), X.numCols() ),// M x N
				Q( X.numRows(), X.numCols() ),// M x N

				Alpha( X.numCols(), X.numCols() ),// N x N
				BetaNeg( X.numCols(), X.numCols() );// N x N. -Beta

	DynamicMatrix<T>	Rs( X.numRows(), X.numCols() ),// M x N
				S( X.numRows(), X.numCols() ),// M x N
				AS( X.numRows(), X.numCols() ),// M x N
				RsAQ( X.numCols(), X.numCols() );// N x N. Rs^t * AQ

	DynamicMatrix<T>	bufNxN( X.numCols(), X.numCols() ),// N x N
				bufMxN( X.numRows(), X.numCols() );// M x N

	QRColumnPivotHouseholder_Solver<T> orth(P, 2);
	//FullPivotLU_Solver<T>	solver(Alpha, 2);
	PivotLU_Solver<T>	solver(Alpha, 2);


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); ++i )
	{
		for( int j=0; j<R.numCols(); ++j )
		{
			T ax = 0;
			for( int k=0; k<A.numCols(); ++k )
				ax += A(i, k) * X(k, j);

			R(i, j) = B(i, j) - ax;	// R0 = B - AX0
			Rs(i, j) = R(i, j);		// R* = R

			P(i, j) = R(i, j);		// P0 = R0
		}
	}


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// Orthogonalize using QR decomposition
		orth.SetCoeff(P);
		for( int i=0; i<Q.numCols(); ++i )
			Q.CopyColumn( orth.Q(), i, i );

		// AQ = A * Qk
		Multiply( AQ, A, Q );

		// Solve (Rs^t * AQ) * Alpha = Rs^t * R
		TransposeMultiply( RsAQ, Rs, AQ );
		solver.SetCoeff( RsAQ );
		TransposeMultiply( bufNxN, Rs, R );// Rs^t * R
		solver.Solve( Alpha, bufNxN );

		// Sk = Rk - AQ * Alpha
		Multiply( bufMxN, AQ, Alpha );
		Subtract( S, R, bufMxN );

		// AS = A * Sk
		Multiply( AS, A, S );

		// z = Tr( AS^t * S ) / Tr( AS^t * AS )
		TransposeMultiply( bufNxN, AS, S );
		T z = Tr( bufNxN );// z = Tr( AS^t * S )
		TransposeMultiply( bufNxN, AS );
		z /= Tr( bufNxN );// z /= Tr( AS^t * AS )

		// X_next += Qk * Alpha + z * Sk
		Multiply( bufMxN, Q, Alpha );// Qk * Alpha
		Add( X, bufMxN );// X_next += Qk * Alpha
		Scale( bufMxN, z, S );//z * Sk
		Add( X, bufMxN );// X_next += z * Sk

		// R_next = Sk - z * ASk
		Scale( bufMxN, z, AS );// z * ASk
		Subtract( R, S, bufMxN );//Sk - z * ASk

		// Check conversion
		TransposeMultiply( bufNxN, R );// R_next^t * R_next
		if( FrobeniusNorm( bufNxN ) < eps )
			break;

		// Solve (Rs^t * AQ) * Beta = -Rs^t * ASk
		TransposeMultiply( bufNxN, Rs, AS );// Rs^t * ASk
		solver.Solve( BetaNeg, bufNxN );// Get NEGATIVE Beta!

		// P_next = R_next + (Q_current - z * AQk) * Beta   ->   P_next = R_next + (-Q_current + z * AQk) * BetaNeg
		AddScaled( bufMxN, Q, (T)-1, AQ, z );// (-Q_current + z * AQk)
		Multiply( P, bufMxN, BetaNeg );// P_next = (-Q_current + z * AQk) * BetaNeg
		Add( P, R );// P_next += R_next

	}// end of k loop

}



// https://arxiv.org/pdf/1104.0737.pdf

template< typename T >
void BlockBiCGSTAB_QR_Precond( const IMatrix<T>& A, IMatrix<T>& X, const IMatrix<T>& B, Solver<T>& precond, int maxiter=100, double eps=1.0e-9 )
{

	DynamicMatrix<T>	R( X.numRows(), X.numCols() ),// M x N
				P( X.numRows(), X.numCols() ),// M x N
				AQ( X.numRows(), X.numCols() ),// M x N
				Q( X.numRows(), X.numCols() ),// M x N

				Alpha( X.numCols(), X.numCols() ),// N x N
				BetaNeg( X.numCols(), X.numCols() );// N x N. -Beta

	DynamicMatrix<T>	Rs( X.numRows(), X.numCols() ),// M x N
				S( X.numRows(), X.numCols() ),// M x N
				AS( X.numRows(), X.numCols() ),// M x N
				RsAQ( X.numCols(), X.numCols() );// N x N. Rs^t * AQ

	DynamicMatrix<T>	bufNxN( X.numCols(), X.numCols() ),// N x N
				bufMxN( X.numRows(), X.numCols() );// M x N

	DynamicMatrix<T>	S_( X.numRows(), X.numCols() );

	QRColumnPivotHouseholder_Solver<T> orth(P, 2);
	//FullPivotLU_Solver<T>	solver(Alpha, 2);
	PivotLU_Solver<T>	solver(Alpha, 2);


	//============== Initialize r and p =============//
	for( int i=0; i<A.numRows(); ++i )
	{
		for( int j=0; j<R.numCols(); ++j )
		{
			T ax = 0;
			for( int k=0; k<A.numCols(); ++k )
				ax += A(i, k) * X(k, j);

			R(i, j) = B(i, j) - ax;	// R0 = B - AX0
			Rs(i, j) = R(i, j);		// R* = R

			P(i, j) = R(i, j);		// P0 = R0
		}
	}


	//================== Iterate ===================//
	for( int k=0; k<maxiter; k++ )
	{
		// Orthogonalize using QR decomposition
		orth.Init(P);
		for( int i=0; i<P.numCols(); ++i )
			P.CopyColumn( orth.Q(), i, i );

		// Preconditioning. Solve M^-1 * Q = P
		precond.Solve( Q, P );

		// AQ = A * Qk
		Multiply( AQ, A, Q );

		// Solve (Rs^t * AQ) * Alpha = Rs^t * R
		TransposeMultiply( RsAQ, Rs, AQ );
		solver.SetCoeff( RsAQ );
		TransposeMultiply( bufNxN, Rs, R );// Rs^t * R
		solver.Solve( Alpha, bufNxN );

		// Sk = Rk - AQ * Alpha
		Multiply( bufMxN, AQ, Alpha );
		Subtract( /*S*/S_, R, bufMxN );

// Solve S = MT  ->  M^-1 * S = T
precond.Solve( S, S_ );


		// AS = A * Sk
		Multiply( AS, A, S );

		// z = Tr( AS^t * S ) / Tr( AS^t * AS )
		TransposeMultiply( bufNxN, AS, S );
		T z = Tr( bufNxN );// z = Tr( AS^t * S )
		TransposeMultiply( bufNxN, AS );
		z /= Tr( bufNxN );// z /= Tr( AS^t * AS )

		// X_next += Qk * Alpha + z * Sk
		Multiply( bufMxN, Q, Alpha );// Qk * Alpha
		Add( X, bufMxN );// X_next += Qk * Alpha
		Scale( bufMxN, z, S );//z * Sk
		Add( X, bufMxN );// X_next += z * Sk

		// R_next = Sk - z * ASk
		Scale( bufMxN, z, AS );// z * ASk
		Subtract( R, S, bufMxN );//Sk - z * ASk

		// Check conversion
		TransposeMultiply( bufNxN, R );// R_next^t * R_next
		if( FrobeniusNorm( bufNxN ) < eps )
			break;

		// Solve (Rs^t * AQ) * Beta = -Rs^t * ASk
		TransposeMultiply( bufNxN, Rs, AS );// Rs^t * ASk
		solver.Solve( BetaNeg, bufNxN );// Get NEGATIVE Beta!

		// P_next = R_next + (Q_current - z * AQk) * Beta   ->   P_next = R_next + (-Q_current + z * AQk) * BetaNeg
		AddScaled( bufMxN, Q, (T)-1, AQ, z );// (-Q_current + z * AQk)
		Multiply( P, bufMxN, BetaNeg );// P_next = (-Q_current + z * AQk) * BetaNeg
		Add( P, R );// P_next += R_next

	}// end of k loop

}







int main()
{
	// Worse condition problem
	double a_[10][10] =
	{
		{ 894,    0,   0,     0,   0,   28,  0,   0,   1000,  70000 },
		{ 0,      5,   13,    5,   0,   0,   0,   0,   0,     0 },
		{ 0,      13,  72,    34,  0,   0,   0,   0,   0,     6500 },
		{ 0,      5,   34,    1,   0,   0,   0,   0,   0,     55 },
		{ 0,      0,   0,     0,   70,  0,   28,  32,  12,    0 },
		{ 28,     0,   0,     0,   0,   87,  20,  0,   33,    0 },
		{ 0,      0,   0,     0,   28,  20,  71,  39,  0,     0 },
		{ 0,      0,   0,     0,   32,  0,   39,  46,  8,     0 },
		{ 1000,   0,   0,     0,   12,  33,  0,   8,   82,    11 },
		{ 70000,  0,   6500,  55,  0,   0,   0,   0,   11,    100 }
	};

	double b_[10][2] =
	{
		//{ 1, 5 },
		//{ 1, 5 },
		//{ 1, 6 },
		//{ 1, 3 },
		//{ 1, 5 },
		//{ 1, 5 },
		//{ 1, 5 },
		//{ 1, 2 },
		//{ 1, -5 },
		//{ 1, 5 }

		//{ 5, 5 },
		//{ 5, 5 },
		//{ 6, 6 },
		//{ 3, 3 },
		//{ 5, 5 },
		//{ 5, 5 },
		//{ 5, 5 },
		//{ 2, 2 },
		//{ -5, -5 },
		//{ 5, 5 }

		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 },
		//{ 1, 1 }

		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 },
		{ 0, 1 }

		//{ 0, 5 },
		//{ 0, 5 },
		//{ 0, 6 },
		//{ 0, 3 },
		//{ 0, 5 },
		//{ 0, 5 },
		//{ 0, 5 },
		//{ 0, 2 },
		//{ 0, -5 },
		//{ 0, 5 }
	};

	//0.0174039
	//1.55702
	//-0.179325
	//-0.890778
	//0.0458074
	//0.13817
	//-0.119308
	//0.139099
	//-0.276422
	//0.00368562

	DynamicMatrix<double>	A( 10, 10, (double *)a_ );
	StaticMatrix<double, 10, 2>	x, b( (double *)b_ );

	tcout << _T("A:\n");
	A.Display();
	tcout << tendl;

	tcout << _T("b:\n");
	b.Display();
	tcout << tendl;

	
	{
		tcout << _T("//========================== BlockCG ==========================//\n");

		Zero(x);

		BlockCG( A, x, b );

		tcout << _T("x:\n");
		x.Display();
		tcout << tendl;

		tcout << _T("Ax - b:\n");
		(A*x - b).Display();
		tcout << tendl;
	}
	
	tcout << tendl;
	
	{
		tcout << _T("//======================= BlockBiCGSTAB =======================//\n");

		Zero(x);

		BlockBiCGSTAB( A, x, b );

		tcout << _T("x:\n");
		x.Display();
		tcout << tendl;

		tcout << _T("Ax - b:\n");
		(A*x - b).Display();
		tcout << tendl;
	}
	
	tcout << tendl;

	{
		tcout << _T("//================== Breakdown-free BlockCG ===================//\n");

		PivotLU_Solver<double> lupivot_precond(A);
		ILU0_Preconditioner<double> ilu0_precond(A);

		Zero(x);

		BreakdownFreeBlockCG( A, x, b, lupivot_precond/*ilu0_precond*/ );

		tcout << _T("x:\n");
		x.Display();
		tcout << tendl;

		tcout << _T("Ax - b:\n");
		(A*x - b).Display();
		tcout << tendl;
	}
	
	tcout << tendl;
	
	{
		tcout << _T("//===================== BlockBiCGSTAB_QR =====================//\n");
	
		Zero(x);

		BlockBiCGSTAB_QR( A, x, b );

		tcout << _T("x:\n");
		x.Display();
		tcout << tendl;

		tcout << _T("Ax - b:\n");
		(A*x - b).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//================= BlockBiCGSTAB_QR_Precond =================//\n");
	
		PivotLU_Solver<double> lupivot_precond(A);
		FullPivotLU_Solver<double> lufullpivot_precond(A);
		ILU0_Preconditioner<double> ilu0_precond(A);
		MILU0_Preconditioner<double> milu0_precond(A);
		ModifiedCholesky_Solver<double> mc_precond(A);
		IncompleteCholesky_Solver<double> ic_precond(A);
		SSORPreconditioner<double> ssor(A,1.0);

		Zero(x);

		BlockBiCGSTAB_QR_Precond( A, x, b, lupivot_precond/*ssor*//*lufullpivot_precond*//*ilu0_precond*//*milu0_precond*//*mc_precond*//*ic_precond*/ );

		tcout << _T("x:\n");
		x.Display();
		tcout << tendl;

		tcout << _T("Ax - b:\n");
		(A*x - b).Display();
		tcout << tendl;
	}

	return 0;
}
