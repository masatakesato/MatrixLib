#ifndef CHOLESKY_H
#define CHOLESKY_H


#include	"Solver.h"
#include	"Substitution.h"



template< typename T >
inline bool Cholesky( IMatrix<T>& L, const IMatrix<T>& A )
{
	if( !A.IsSquare() || !A.IsSameSize(L) )
		return false;

	for( int i=0; i<L.numRows(); ++i )
	{
		//============ Fill non-diagonal elements =============//
		for( int j=0; j<i; ++j )
		{
			T s = A(i, j);
			for( int k=0; k<j; ++k )	s -= L(i, k) * L(j, k);

			L(i, j) = s / L(j, j);// Calculate L[i][j]
		}// end of j loop

		//============== Fill diagonal elements ==============//
		T s = A(i, i);
		for( int k=0; k<i; ++k )	s -= pow( L(i, k), 2 );

		L(i, i) = sqrt(s); // Calculate L[i][i]

	}// end of i loop

	return true;
}



// http://fornext1119.hatenablog.com/entry/2014/07/30/131142

template< typename T >
inline bool ModifiedCholesky( IMatrix<T>& L, IMatrix<T>& D, const IMatrix<T>& A )
{
	if( !(L.IsSameSize(D) && L.IsSameSize(A) && L.IsSquare()) )
		return false;

	D(0, 0) = A(0, 0);
	L(0, 0) = (T)1;

	for( int i=1; i<L.numRows(); ++i )
	{
		for( int j=0; j<i; ++j )
		{
			T l_ij = A(i, j);
			for( int k=0; k<j; ++k )
				l_ij -= L(i, k) * D(k, k) * L(j, k);

			L(i, j) = l_ij / D(j, j);
		}

		T d_ii = A(i, i);
		for( int k=0; k<i; ++k )
			d_ii -= L(i, k) * L(i, k) * D(k, k);

		D(i, i) = d_ii;
		L(i, i) = (T)1;

	}// end of i loop

	return true;
}



template< typename T >
inline bool IncompleteCholesky( IMatrix<T>& L, IMatrix<T>& D, const IMatrix<T>& A )
{
	if( !(L.IsSameSize(D) && L.IsSameSize(A) && L.IsSquare()) )
		return false;

	D(0, 0) = A(0, 0);
	L(0, 0) = (T)1;

	for( int i=1; i<L.numRows(); ++i )
	{
		for( int j=0; j<i; ++j )
		{
			if( fabs( A(i, j) ) < 1.0e-10 )	continue;

			T l_ij = A(i, j);
			for( int k=0; k<j; ++k )
				l_ij -= L(i, k) * D(k, k) * L(j, k);

			L(i, j) = l_ij / D(j, j);
		}

		T d_ii = A(i, i);
		for( int k=0; k<i; ++k )
			d_ii -= L(i, k) * L(i, k) * D(k, k);

		D(i, i) = d_ii;
		L(i, i) = (T)1;

	}// end of i loop

	return true;
}



// Inverse_Cholesky
// https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
template< typename T >
bool Inverse_Cholesky( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( !(in.IsSameSize(out) && in.IsSquare()) )
		return false;

	DynamicMatrix<T> L_inv( in.numRows(), in.numCols() );//, L_inv_t( in.numRows(), in.numCols() );

	// Use out as decomposed_cholesky temporarily.
	Cholesky( out, in );

	if( out.IsFinite()==false )// Abort if A is NOT positive definite matrix
	{
		tcerr << "Cholesky_Solver::Unable to decompose... Matrix is not positive definite.\n";
		return false;
	}

	// Calculate L_inv and U_inv
	for( int n=0; n<L_inv.numCols(); ++n )
	{
		// Solve L * L_inv = I with forward substitution
		for( int i=n; i<out.numRows(); ++i )
		{
			L_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
			for( int j=0; j<i; ++j )
				L_inv(i, n) -= out(i, j) * L_inv(j, n);
			L_inv(i, n) /= out(i, i);
		}
	}

	//out.Display();
	//L_inv.Display();
	//L_inv_t.Display();

	// in^-1 = L_inv_t * L_inv
	TransposeMultiply( out, L_inv );

	return true;
}




template< typename T >
inline bool IsPositiveDefinite( const IMatrix<T>& A )
{
	Matrix<T> L( A.numRows(), A.numCols() );
	Cholesky( L, A );

	return L.IsFinite();
}





template< typename T >
class Cholesky_Solver : public Solver<T>
{
public:

	Cholesky_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~Cholesky_Solver()
	{
		m_L.Release();
		m_y.Release();

		L_inv.Release();

		m_Inverse.Release();
		m_bInvUpdated = false;
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		m_L.Init( A.numRows(), A.numCols() );

		if( SetCoeff( A ) ==false )
		{
			m_L.Release();
			return;
		}

		m_y.Init( A.numCols(), batchSize );
		L_inv.Init( A.numRows(), A.numCols() );
		m_Inverse.Init( A.numRows(), A.numCols() );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		m_bInvUpdated = false;

		bool result = Cholesky( m_L, A );

		if( m_L.IsFinite()==false )// Abort if A is NOT positive definite matrix
		{
			tcout << "Cholesky_Solver::Unable to decompose... Matrix is not positive definite.\n";
			return false;
		}

		return result;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_L.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_L.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_L.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve L * y = b using forward substitution( Equivalent to SUbstitution::ForwardSubstitution( m_L, m_y, b ) )
			for( int i=0; i<m_L.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; j++ )
					m_y(i, tid) -= m_L(i, j) * m_y(j, tid);
				m_y(i, tid) /= m_L(i, i);
			}

			// Solve Lt * x = y using backward substitution
			for( int i=m_L.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_L.numCols(); ++j )
					x(i, k) -= m_L(j, i) * x(j, k);// access m_L using transosed index(j, i)
				x(i, k) /= m_L(i, i);
			}
		}

	}


	const IMatrix<T>& Inverse()
	{
		if(m_bInvUpdated==true)
			return m_Inverse;

		// Calculate L_inv
		for( int n=0; n<L_inv.numCols(); ++n )
		{
			// Solve L * L_inv = I with forward substitution
			for( int i=n; i<m_L.numRows(); ++i )
			{
				L_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=0; j<i; ++j )
					L_inv(i, n) -= m_L(i, j) * L_inv(j, n);
				L_inv(i, n) /= m_L(i, i);
			}
		}

		// in^-1 = L_inv_t * L_inv
		TransposeMultiply( m_Inverse, L_inv );

		m_bInvUpdated = true;

		return m_Inverse;
	}



private:

	DynamicMatrix<T>	m_L;	// Composed L matrix
	DynamicMatrix<T>	m_y;	// Intermediate variables.

	// temporal buffer for inverse matrix calculation
	DynamicMatrix<T>	L_inv;

	// Inverse Matrix of A
	DynamicMatrix<T>	m_Inverse;
	bool		m_bInvUpdated;
};




template< typename T >
class ModifiedCholesky_Solver : public Solver<T>
{
public:

	ModifiedCholesky_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~ModifiedCholesky_Solver()
	{
		m_LD.Release();
		m_Lt.Release();
		m_y.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		m_LD.Init( A.numRows(), A.numCols() );
		m_Lt.Init( A.numRows(), A.numCols() );

		if( SetCoeff( A )==false )
		{
			m_LD.Release();
			m_Lt.Release();
			return;
		}

		m_y.Init( A.numCols(), batchSize );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		DynamicMatrix<T> l( A.numRows(), A.numCols() ), d( A.numRows(), A.numCols() );
		bool result = ModifiedCholesky( l, d, A );

		if( l.IsFinite()==false )// Abort if A is NOT positive definite matrix
		{
			tcout << "ModifiedCholesky_Solver::Unable to decompose... A is not positive definite matrix.\n";
			return false;
		}

		Multiply( m_LD, l, d );
		Transpose( m_Lt, l );

		return result;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_LD.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_LD.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_LD.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve (L*D) * y = b 
			//ForwardSubstitution( m_LD, m_y, b );
			for( int i=0; i<m_LD.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_LD(i, j) * m_y(j, tid);
				m_y(i, tid) /= m_LD(i, i);
			}

			// Solve Lt * x = y
			//BackwardSubstitution( m_Lt, x, m_y );
			for( int i=m_Lt.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Lt.numCols(); ++j )
					x(i, k) -= m_Lt(i, j) * x(j, k);
				x(i, k) /= m_Lt(i, i);
			}
		}

	}



private:

	DynamicMatrix<T>	m_LD;	// L*D matrix
	DynamicMatrix<T>	m_Lt;	// Lt matrix
	DynamicMatrix<T>	m_y;	// Intermediate variables.

};




template< typename T >
class IncompleteCholesky_Solver : public Solver<T>
{
public:

	IncompleteCholesky_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~IncompleteCholesky_Solver()
	{
		m_LD.Release();
		m_Lt.Release();
		m_y.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		m_LD.Init( A.numRows(), A.numCols() );
		m_Lt.Init( A.numRows(), A.numCols() );

		if( SetCoeff( A )==false )
		{
			m_LD.Release();
			m_Lt.Release();
			return;
		}	

		m_y.Init( A.numCols(), batchSize );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		DynamicMatrix<T> l( A.numRows(), A.numCols() ), d( A.numRows(), A.numCols() );
		bool result = IncompleteCholesky( l, d, A );

		if( l.IsFinite()==false )// Abort if A is NOT positive definite matrix
		{
			tcout << "IncompleteCholesky_Solver::Unable to decompose... A is not positive definite matrix.\n";
			return false;
		}

		Multiply( m_LD, l, d );
		Transpose( m_Lt, l );

		return result;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_LD.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_LD.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_LD.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve (L*D) * y = b 
			//ForwardSubstitution( m_LD, m_y, b );
			for( int i=0; i<m_LD.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_LD(i, j) * m_y(j, tid);
				m_y(i, tid) /= m_LD(i, i);
			}

			// Solve Lt * x = y
			//BackwardSubstitution( m_Lt, x, m_y );
			for( int i=m_Lt.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Lt.numCols(); ++j )
					x(i, k) -= m_Lt(i, j) * x(j, k);
				x(i, k) /= m_Lt(i, i);
			}
		}

	}



private:

	DynamicMatrix<T>	m_LD;	// L*D matrix
	DynamicMatrix<T>	m_Lt;	// Lt matrix
	DynamicMatrix<T>	m_y;	// Intermediate variables.

};





//##############################################################################//
//							Experimental implementation							//
//##############################################################################//

/*

// Solve Ax = b by Cholesky decomposition.
template< typename T >
inline void SolveCholesky( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b )
{
	// L * Lt * x = b
	Matrix<T> L( A.numRows(), A.numCols() ), Lt(  A.numRows(), A.numCols() );
	Cholesky( L, A );
	Transpose( Lt, L );

	// Lt * x = y

	// Solve L * y = b
	Matrix<T> y( x.numRows(), x.numCols() );
	ForwardSubstitution( L, y, b );

	// Solve Lt * x = y
	BackwardSubstitution( Lt, x, y );

}



// Solve Ax = b by Modified Cholesky decomposition
template< typename T >
inline void SolveModifiedCholesky( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b )
{
	// L * Lt * x = b
	Matrix<T> L( A.numRows(), A.numCols() ), D( A.numRows(), A.numCols() ), Lt( A.numRows(), A.numCols() );
	ModifiedCholesky( L, D, A );
	Transpose( Lt, L );

	// Lt * x = y

	// Solve L * D * y = b
	Matrix<T> y( x.numRows(), x.numCols() );
	ForwardSubstitution( L*D, y, b );

	// Solve Lt * x = y
	BackwardSubstitution( Lt, x, y );

}




// Solve Ax = b by Modified Incomplete Cholesky decomposition
template< typename T >
inline void SolveIncompleteCholesky( const Matrix<T>& A, Matrix<T>& x, const Matrix<T>& b )
{
	// L * Lt * x = b
	Matrix<T> L( A.numRows(), A.numCols() ), D( A.numRows(), A.numCols() ), LD( A.numRows(), A.numCols() ), Lt( A.numRows(), A.numCols() );
	IncompleteCholesky( L, D, A );
	Transpose( Lt, L );
	Multiply( LD, L, D );

	// Lt * x = y

	// Solve L * D * y = b
	Matrix<T> y( x.numRows(), x.numCols() );
	ForwardSubstitution( LD, y, b );

	// Solve Lt * x = y
	BackwardSubstitution( Lt, x, y );

}

*/


#endif /* CHOLESKY_H */