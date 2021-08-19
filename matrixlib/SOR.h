#ifndef SOR_H
#define	SOR_H


#include	"Solver.h"
#include	"Substitution.h"



template< typename T >
inline void SOR( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, double omega, int maxiter=100, double eps=1.0e-9 )
{

	for( int k=0; k<maxiter; ++k )
	{
		T err = 0;

		for( int i=0; i<x.numRows(); ++i )
		{
			// Calculate x_k+1
			T tmp = x(i, 0);
			T x_i = b(i, 0);

			for( int j=0; j<A.numCols(); ++j )
			{
				if( i==j ) continue;
				x_i -= A(i, j) * x(j, 0);
			}
			x_i /= A(i, i);
			x(i, 0) += (T)omega * ( x_i - x(i, 0) );

			err += abs( tmp - x_i );

			//tcout << x_i << tendl;
		}// end of i loop


		//x.Display();
		 // Check convergence
		if( err < eps )	break;

	}// end of k loop

}





template< typename T >
inline void SSOR( const IMatrix<T>& A, IMatrix<T>& x, const IMatrix<T>& b, double omega, int maxiter=100, double eps=1.0e-9 )
{

	for( int k=0; k<maxiter; ++k )
	{
		T err = 0;

		// Forward
		for( int i=0; i<x.numRows(); ++i )
		{
			// Calculate x_k+1
			T tmp = x(i, 0);
			T x_i = b( i, 0 );

			for( int j=0; j<A.numCols(); ++j )
			{
				if( i==j ) continue;
				x_i -= A(i, j) * x(j, 0);
			}
			x_i /= A(i, i);
			x(i, 0) += (T)omega * ( x_i - x(i, 0) );

			err += abs( tmp - x_i );

		}// end of i loop

		// Backward
		for( int i=x.numRows()-1; i>=0; --i )
		{
			// Calculate x_k+1
			T tmp = x(i, 0);
			T x_i = b( i, 0 );

			for( int j=0; j<A.numCols(); ++j )
			{
				if( i==j ) continue;
				x_i -= A(i, j) * x(j, 0);
			}
			x_i /= A(i, i);
			x(i, 0) += (T)omega * ( x_i - x(i, 0) );

			err += abs( tmp - x_i );

		}// end of i loop



		 // Check convergence
		if( err < eps )	break;

	}// end of k loop

}




template< typename T >
class SOR_Solver : public Solver<T>
{
public:

	SOR_Solver( const IMatrix<T>& A, double omega, int batchSize=1 )
	{
		Init( A, batchSize );

		maxiter	= 100;
		this->omega	= omega;
		eps		= 1.0e-9;

	}


	~SOR_Solver()
	{
		m_A.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_A.Init( m, n );
		this->batchSize = batchSize;

		SetCoeff( A );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		if( !m_A.IsSameSize(A) )
			return false;

		m_A.Copy( A );

		return true;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_A.IsEmpty() )
		{
			this->batchSize = batchSize;
			return true;
		}
		return false;
	}



	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_A.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			MatrixView<T> x_( x, 0, k, b.numRows(), 1 );
			MatrixView<T> b_( b, 0, k, b.numRows(), 1 );

			SOR( m_A, x_, b_, omega, maxiter, eps );
		}
	}




	#ifdef _OPENMP
	
	void Solve_Parallel( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		#pragma omp parallel for shared(x) num_threads(batchSize)
		for( int k=0; k<b.numCols(); ++k )
		{
			MatrixView<T> x_( x, 0, k, b.numRows(), 1 );
			MatrixView<T> b_( b, 0, k, b.numRows(), 1 );

			SOR( m_A, x_, b_, omega, maxiter, eps );
		}
	}
	
	#endif



private:

	DynamicMatrix<T>	m_A;

	int		maxiter;
	double	omega;
	double	eps;
	int		batchSize;

};








// http://www.numericalmethod.com/javadoc/suanshu/com/numericalmethod/suanshu/algebra/linear/matrix/doubles/matrixtype/sparse/solver/iterative/preconditioner/SSORPreconditioner.html
template< typename T >
class SSORPreconditioner : public Solver<T>
{
public:

	SSORPreconditioner( const IMatrix<T>& A, const T& omega, int batchSize=1 )
	{
		maxiter	= 100;
		eps		= 1.0e-9;

		Init( A, omega, batchSize );
	}


	~SSORPreconditioner()
	{
		
	}


	void Init( const IMatrix<T>& A, const T& omega, int batchSize=1 )
	{
		if( A.IsSymmetric() )// TODO: Check if all diagonal elements are NON-ZERO
		{
			int m = A.numRows();
			int n = A.numCols();

			m_y.Init( n, batchSize );


			m_Dinv.Init( m, n );

			m_Lower.Init( m, n );
			m_Upper.Init( n, m );


			L.Init( m, n );
			D.Init( m, n );


			SetCoeff( A, omega );
		}
	}


	bool SetCoeff( const IMatrix<T>& A, const T& omega=1.06 )
	{
		if( !m_Lower.IsSameSize(A) )
			return false;

		int m = A.numRows();
		int n = A.numCols();
		
		//============== On omega version =================//
		/*
		for( int i=0; i<m; ++i )
		{
			D(i, i) = A(i, i);
			m_Dinv(i, i) = 1 / A(i, i);// D_inv

			for( int j=0; j<i; ++j )
				L(i, j) = A(i, j);
		}
		*/

		//// M = (D + L) * D^inv * (D + L)^t
		//for( int i=0; i<m; ++i )
		//{
		//	// ( D + L ) * D_inv
		//	for( int j=0; j<=i; ++j )
		//		m_Lower(i, j) = A(i, j) / A(j, j);

		//	// ( D + L )^t
		//	for( int j=0; j<=i; ++j )
		//		m_Upper(j, i) = A(i, j);
		//}
		//
		//(m_Lower * m_Upper).Display();


		//================== omega version =====================//
		//for( int i=0; i<m; ++i )
		//{
		//	D(i, i) = A(i, i) / (T)omega;
		//	m_Dinv(i, i) = (T)omega / A(i, i);// D_inv

		//	for( int j=0; j<i; ++j )
		//		L(i, j) = A(i, j);
		//}

		//auto DLt = D + L;
		//DLt.Transpose();
		//auto M = (D + L) * m_Dinv * DLt;
		//Scale(M, 1/(2-(T)omega));
		//M.Display();


		T coeff = T( omega / (2-omega) );

		for( int i=0; i<m; ++i )
		{
			for( int j=0; j<=i; ++j )
			{
				// ( D/omega + L )^t
				m_Upper(j, i) = A(i, j) / T(i==j ? omega : 1);

				// ( D/omega + L ) * 1/(2-omega) * (D/omega)^-1
				m_Lower(i, j) = m_Upper(j, i) * coeff / A(j, j);
			}
		}
		
		//(m_Lower * m_Upper).Display();
		m_Lower.Display();
		m_Upper.Display();

		return true;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Lower.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Lower.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Lower.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve L * y = b using forward substitution( Equivalent to SUbstitution::ForwardSubstitution( m_L, m_y, b ) )
			for( int i=0; i<m_Lower.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Lower(i, j) * m_y(j, tid);
				m_y(i, tid) /= m_Lower(i, i);
			}

			// Solve Lt * x = y using backward substitution
			for( int i=m_Upper.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Upper.numCols(); ++j )
					x(i, k) -= m_Upper(i, j) * x(j, k);
				x(i, k) /= m_Upper(i, i);
			}
		}

	}


	#ifdef _OPENMP
	
	void Solve_Parallel( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		#pragma omp parallel for shared(x) num_threads(m_y.numCols())
		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve L * y = b using forward substitution( Equivalent to SUbstitution::ForwardSubstitution( m_L, m_y, b ) )
			for( int i=0; i<m_Lower.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Lower(i, j) * m_y(j, tid);
				m_y(i, tid) /= m_Lower(i, i);
			}

			// Solve Lt * x = y using backward substitution
			for( int i=m_Upper.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Upper.numCols(); ++j )
					x(i, k) -= m_Upper(i, j) * x(j, k);
				x(i, k) /= m_Upper(i, i);
			}
		}
	}
	
	#endif



private:

	DynamicMatrix<T>	L, D, m_Dinv;
	DynamicMatrix<T>	m_Lower,	// (L + D) * D^inv
				m_Upper;	// (L + D) * t

	DynamicMatrix<T>	m_y;	// Intermediate variables.

	int		maxiter;
	double	eps;


	using Solver<T>::Init;
	bool SetCoeff( const IMatrix<T>& A ){ return false ;}
};




#endif // !SOR_H
