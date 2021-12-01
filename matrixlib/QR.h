#ifndef QR_H
#define	QR_H

#include	<oreore/mathlib/MathLib.h>
#include	<oreore/algorithm/Algorithm.h>

//#include	"Matrix.h"
#include	"Substitution.h"
#include	"Solver.h"



// QR decomposition using modified gram schmidt orthogonalization
template< typename T >
inline bool QR_GramSchmidt( IMatrix<T>& v, IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A )
{
	int m = A.numRows();// vector dimension
	int n = A.numCols();// num of vectors

	// vj = xj
	v.Copy(A);

	for( int j=0; j<n; ++j )
	{
		// rjj = ||vj||2
		T& rjj = R(j, j);
		rjj = 0;
		for( int l=0; l<m; ++l )
			rjj += pow( v(l, j), 2 );
		rjj = sqrt(rjj);

		// qj = vj / rjj
		for( int l=0; l<m; ++l )
			Q(l, j) = v(l, j) / rjj;

		for( int k=j+1; k<n; ++k )
		{
			// rjk = (qj,  vk)
			T& rjk = R(j, k);
			rjk = 0;
			for( int l=0; l<m; ++l )
				rjk += Q(l, j) * v(l, k);

			// vk = vk - rjk * qj
			for( int l=0; l<m; ++l )
				v(l, k) -= rjk * Q(l, j);
		}// end of k loop

	}// end of i loop

	return true;
}



template< typename T >
inline bool QR_GramSchmidt( IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A )
{
	DynamicMatrix<T> v(A);
	return QR_GramSchmidt( v, Q, R, A );
}





// https://media.cheggcdn.com/media/89b/89becc04-e419-4219-9c1d-b27381cad344/phpy83VTr.png
// https://media.cheggcdn.com/media/a6c/a6c9e6df-b814-43fd-af61-29c722190598/phpSJ69l6.png	<-
// http://zakii.la.coocan.jp/signal/38_qr_householder.htm

// QR decomposition using Householder transformation
template< typename T >
inline bool QR_Householder( IMatrix<T>& v, IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A )
{
	int m = A.numRows();// vector dimension
	int n = A.numCols();// num of vectors

	R.Copy(A);
	Identity(Q);

	for( int k=0; k<n; ++k )
	{
		T vk_len = 0;
		// 6: copy R( k:m-1, k ) to v( k:m-1, k )
		for( int l=k; l<m; ++l )
		{
			v(l, 0) = R(l, k);
			vk_len += pow( v(l, 0), 2 ) * T(l>k);// accumulate squared distance except diagonal element.
		}

		// 7:
		v(k, 0) = v(k, 0) + Sign( v(k, 0) ) * sqrt( pow(v(k, 0), 2) + vk_len );

		// 8: normalize vk
		vk_len = 1 / sqrt( pow(v(k, 0), 2) + vk_len );
		for( int l=k; l<m; ++l )	v(l, 0) *= vk_len;

		// 9:	R(k:m-1, k:n-1) = R(k:m-1, k:n-1) - 2 * vk * ( vk^t * R(k:m-1, k:n-1) )
		for( int i=k; i<n; ++i )
		{
			// multiply v's column k by R's column i first.
			T vkr = 0;
			for( int j=k; j<m; ++j )
				vkr += v(j, 0) * R(j, i);

			// update R's column k using vkr.
			for( int j=k; j<m; ++j )
				R(j, i) = R(j, i) - 2 * v(j, 0) * vkr;
		}

		// 10:  Q(k, k:m-1) = Q(k, k:m-1) - 2 * vk * ( vk^t * Q(k, k:m-1) )
		for( int i=0; i<m; ++i )
		{
			// multiply v's column k by Q's column i first.
			T vkq = 0;
			for( int j=k; j<m; ++j )
				vkq += v(j, 0) * Q(i, j);//Q(j, i);//
				// Fetch Q's transposed element value..........(1)

			// update Q's column i using vkq.
			for( int j=k; j<m; ++j )
				Q(i, j) -= 2 * v(j, 0) * vkq;//Q(j, i) -= 2 * v(j, 0) * vkq;//
				// Update Q's transposed element...............(2)
		}

	}// end of k loop

	//Q.Transpose();// Post transpose can be avoided by (1) and (2)

	return true;
}



template< typename T >
inline bool QR_Householder( IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A )
{
	DynamicMatrix<T> v(A);
	return QR_Householder( v, Q, R, A );
}




// https://www.arxiv-vanity.com/papers/1804.05138/
// http://faculty.nps.edu/borges/Teaching/MA3046/Matlab/qrlsdiary.html

template< typename T >
inline int QR_ColumnPivotHouseholder( IMatrix<T>& v, OreOreLib::Array<T>& rs, /*DynamicMatrix<T>& P*/OreOreLib::Array<int>& P, IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A, double epsilon=1.0e-8 )
{
	int k = 0;
	int m = A.numRows();// vector dimension
	int n = A.numCols();// num of vectors

	//OreOreLib::Array<T> rs(n);// l2-norms of A's each column

	for( int i=0; i<n; i++ ) P[i] = i;//Identity(P);
	Identity(Q);
	R.Copy(A);

	// Initialize rs[1:n]
	for( int s=0; s<n; ++s )
		rs[s] = DotCC( R, R, s, s );


	for( k=0; k<n; ++k )
	{
		//================= Column pivoting ==================//
		//for( int i=0; i<n; ++i )
		//	tcout << rs[i] << ", ";
		//tcout << tendl;

		// find j = argmax from rs[k:n]
		int j = (int)OreOreLib::ArgMax( rs.begin()+k, rs.end() ) + k;
		if( rs[j] < epsilon )
			break;
		// Swap row k and j of rs
		rs.Swap( k, j );

		// Swap col k and j of A
		R.SwapColumn( k, j );

		// Update permutation by swapping col k and j
		//P.SwapColumn( k, j );
		P.Swap( k, j );

		//============= Householder reflection ===============//

		// 6: copy R( k:m-1, k ) to v( k:m-1, k )
		v.CopyColumn( R, k, 0 );

		T vk_len = rs[k] - pow( v(k, 0), 2 );// accumulate squared distance except diagonal element.

											 // 7:
		v(k, 0) = v(k, 0) + Sign( v(k, 0) ) * sqrt( pow(v(k, 0), 2) + vk_len );

		// 8: normalize vk
		T vk_len_inv = 1 / sqrt( pow(v(k, 0), 2) + vk_len );
		for( int l=k; l<m; ++l )
			v(l, 0) *= vk_len_inv;

		// 9:	R(k:m-1, k:n-1) = R(k:m-1, k:n-1) - 2 * vk * ( vk^t * R(k:m-1, k:n-1) )
		for( int i=k; i<n; ++i )
		{
			// multiply v's column k by R's column i first.
			T vkr = 0;
			for( int j=k; j<m; ++j )
				vkr += v(j, 0) * R(j, i);

			// update R's column k using vkr.
			for( int j=k; j<m; ++j )
				R(j, i) = R(j, i) - 2 * v(j, 0) * vkr;
		}

		// 10:  Q(k, k:m-1) = Q(k, k:m-1) - 2 * vk * ( vk^t * Q(k, k:m-1) )
		for( int i=0; i<m; ++i )
		{
			// multiply v's column k by Q's column i first.
			T vkq = 0;
			for( int j=k; j<m; ++j )
				vkq += v(j, 0) * Q(i, j);//Q(j, i);//
				// Fetch Q's transposed element value..........(1)

			// update Q's column i using vkq.
			for( int j=k; j<m; ++j )
				Q(i, j) -= 2 * v(j, 0) * vkq;//Q(j, i) -= 2 * v(j, 0) * vkq;//
				// Update Q's transposed element...............(2)
		}

		//================= Column pivoting ==================//

		// Update rs[k+1:n] for next loop
		for( int s=k+1; s<n; ++s )
			rs[s] -= pow( R(k, s), 2 );

	}// end of k loop

	//Q.Transpose();// Post transpose can be avoided by (1) and (2)
	
	return k;
}



template< typename T >
inline int QR_ColumnPivotHouseholder( OreOreLib::Array<int>& P, IMatrix<T>& Q, IMatrix<T>& R, const IMatrix<T>& A, double epsilon=1.0e-8 )
{
	DynamicMatrix<T> v( A.numRows(), 1 );
	OreOreLib::Array<T> rs( A.numCols() );
	return QR_ColumnPivotHouseholder( v, rs, P, Q, R, A, epsilon );
}





template< typename T >
class QRGramSchmidt_Solver : public Solver<T>
{
public:

	QRGramSchmidt_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~QRGramSchmidt_Solver()
	{
		m_y.Release();

		m_v.Release();
		m_Q.Release();
		m_R.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_y.Init( n, batchSize );

		m_v.Init( m, n );
		m_Q.Init( m, m );
		m_R.Init( m, n );

		SetCoeff(A);
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		return QR_GramSchmidt( m_v, m_Q, m_R, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_R.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_R.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_R.numRows()) )
			return;

		// Solve Ax=b as (Q*R)x = b
		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();
			
			// Solve Q*y=b(y=R*x) -> y=Q^t*b
			for( int i=0; i<m_y.numRows(); ++i )
			{
				m_y(i, tid)	= 0;
				for( int j=0; j<m_Q.numRows(); ++j )
					m_y(i, tid) += m_Q(j, i) * b(j, k);// Multiply m_Q^t by k-th column of b
			}

			//BackwardSubstitution( m_R, x, m_y );// Solve y=R*x
			for( int i=m_R.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_R.numCols(); ++j )
					x(i, k) -= m_R(i, j) * x(j, k);
				x(i, k) /= m_R(i, i);
			}
		}
	}


	const DynamicMatrix<T>& Q() const
	{
		return m_Q;
	}


	const DynamicMatrix<T>& R() const
	{
		return m_R;
	}



private:

	DynamicMatrix<T>	m_y;	// Intermediate variables.

	DynamicMatrix<T>	m_v, m_Q, m_R;

};



template< typename T >
class QRHouseholder_Solver : public Solver<T>
{
public:

	QRHouseholder_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~QRHouseholder_Solver()
	{
		m_y.Release();

		m_v.Release();
		m_Q.Release();
		m_R.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_y.Init( n, batchSize );

		m_v.Init( m, 1 );
		m_Q.Init( m, m );
		m_R.Init( m, n );

		SetCoeff(A);
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		return QR_Householder( m_v, m_Q, m_R, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_R.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_R.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_R.numRows()) )
			return;

		// Solve Ax=b as (Q*R)x = b
		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve Q*y=b(y=R*x) -> y=Q^t*b
			for( int i=0; i<m_y.numRows(); ++i )
			{
				m_y(i, tid)	= 0;
				for( int j=0; j<m_Q.numRows(); ++j )
					m_y(i, tid) += m_Q(j, i) * b(j, k);// Multiply m_Q^t by k-th column of b
			}

			//BackwardSubstitution( m_R, x, m_y );// Solve y=R*x
			for( int i=m_R.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_R.numCols(); ++j )
					x(i, k) -= m_R(i, j) * x(j, k);
				x(i, k) /= m_R(i, i);
			}
		}
	}


	const DynamicMatrix<T>& Q() const
	{
		return m_Q;
	}


	const DynamicMatrix<T>& R() const
	{
		return m_R;
	}



private:

	DynamicMatrix<T>	m_y;	// Intermediate variables.
	DynamicMatrix<T>	m_v, m_Q, m_R;

};




template< typename T >
class QRColumnPivotHouseholder_Solver : public Solver<T>
{
public:

	QRColumnPivotHouseholder_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~QRColumnPivotHouseholder_Solver()
	{
		m_Rank = 0;

		m_y.Release();
		m_z.Release();

		m_v.Release();
		m_Rs.Release();
		m_Pivot.Release();//m_P.Release();
		m_Q.Release();
		m_R.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_y.Init( n, batchSize );
		m_z.Init( n, batchSize );

		m_v.Init( m, 1 );
		m_Rs.Init( n );
		m_Pivot.Init( n );//m_P.Init( n, n );
		m_Q.Init( m, m );
		m_R.Init( m, n );

		SetCoeff(A);
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		m_Rank = QR_ColumnPivotHouseholder( m_v, m_Rs, /*m_P*/m_Pivot, m_Q, m_R, A );
		return m_Rank > 0;
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_R.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_R.numCols(), batchSize );
			return true;
		}
		return false;
	}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_R.numRows()) )
			return;

		// Solve Ax=b as (Q*R)x = b
		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Solve Q*y=b(y=R*x) -> y=Q^t*b
			for( int i=0; i<m_y.numRows(); ++i )
			{
				m_y(i, tid)	= 0;
				for( int j=0; j<m_Q.numRows(); ++j )
					m_y(i, tid) += m_Q(j, i) * b(j, k);// Multiply m_Q^t by k-th column of b
			}

			//BackwardSubstitution( m_R, z, m_y );// Solve y=R*z
			for( int i=m_R.numRows()-1; i>=0; --i )
			{
				m_z(i, tid) = m_y(i, tid);
				for( int j=i+1; j<m_R.numCols(); ++j )
					m_z(i, tid) -= m_R(i, j) * m_z(j, tid);
				m_z(i, tid) /= m_R(i, i);

				// x = Perm(z)
				x(m_Pivot[i], k) = m_z(i, tid);
			}

		}

	}


	const int Rank() const
	{
		return m_Rank;
	}


	//const DynamicMatrix<T>& P() const
	//{
	//	return m_P;
	//}


	const DynamicMatrix<T>& Q() const
	{
		return m_Q;
	}


	const DynamicMatrix<T>& R() const
	{
		return m_R;
	}

	const OreOreLib::Array<int>& Pivot() const
	{
		return m_Pivot;
	}



private:

	DynamicMatrix<T>	m_y, m_z;	// Intermediate variables.

	int			m_Rank;
	DynamicMatrix<T>	m_v, /*m_P,*/ m_Q, m_R;

	OreOreLib::Array<T>		m_Rs;// l2-norms of A's each column
	OreOreLib::Array<int>	m_Pivot;
};




#endif // !QR_H