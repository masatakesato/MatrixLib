#ifndef LU_H
#define	LU_H



//#include	"Matrix.h"
#include	"Solver.h"
//#include	"Substitution.h"


//L = {
//	{ 1,      0,      0, ... }
//	{ l(1,0), 1,      0, ... }
//	{ l(2,0), l(2,1), 1, ... }
//	{ ...
//}


//U = {
//	{ u(0,0), u(0,1), u(0,2), ... }
//	{ 0,      u(1,1), u(1,2), ... }
//	{ 0,      0,      u(2,2), ... }
//	{ ...
//} 



// http://nkl.cc.u-tokyo.ac.jp/11e/05-Solver/solver.pdf

template< typename T >
inline bool LU( IMatrix<T>& a )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();

	for( int k=0; k<m-1; ++k )
	{
		for( int i=k+1; i<m; ++i )
		{
			if( abs(a(k, k)) < 1.0e-10 )
				continue;

			a(i, k) = a(i, k) / a(k, k);
			for( int j=k+1; j<n; ++j )
				a(i, j) -= a(i, k) * a(k, j);
		}
	}// end of k loop

	return true;
}



template< typename T >
inline bool LU( IMatrix<T>&lu, const IMatrix<T>& a )
{
	return lu.Copy(a) ? LU(lu) : false;
}



// http://www.ed.u-tokai.ac.jp/takakura/jp/takakura/book_SuuchiKeisan/RC04.pdf
template< typename T >
inline bool LU_Pivot( IMatrix<T>& a, OreOreLib::Array<int>& p )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();
	int r;
	T amax;

	// Initialize pivot indices
	for( int i=0; i<m; ++i )
		p[i] = i;

	// Forward Substitution with Pivoting
	for( int k=0; k<m-1; ++k )
	{
		//================ Pivotting ================//
		amax = abs( a(k, k) );
		r = k;
		for( int i=k+1; i<m; ++i ) // search r such that a_{r,k}=max a_{i,k} i=k,..n
		{
			if( abs(a(i, k)) > amax )
			{
				amax=abs(a(i, k));
				r=i;
			}
		}

		if( r!=k ) // exchange column r and column k
		{
			p.Swap( r, k );
			a.SwapRow( r, k );
		}

		if( a(k, k)==0 )
			continue;
				
		//============== LU decomposition ============//
		for( int i=k+1; i<m; ++i )
		{
			a(i, k) = a(i, k) / a(k, k);

			for( int j=k+1; j<n; ++j )
				a(i, j) -= a(i, k) * a(k, j);
		}

	}// end of k loop

	return true;
}



template< typename T >
inline bool LU_Pivot( IMatrix<T>&lu, OreOreLib::Array<int>& p, const IMatrix<T>& a )
{
	return lu.Copy(a) ? LU_Pivot(lu, p) : false;
}







template< typename T >
inline bool LU_FullPivot( IMatrix<T>& a, OreOreLib::Array<int>& p, OreOreLib::Array<int>& q )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();
	int mu, lambda;
	
	T amax;

	// Initialize pivot indices
	for( int i=0; i<m; ++i )
	{
		p[i] = i;
		q[i] = i;
	}

	// Forward Substitution with Pivoting
	for( int k=0; k<m-1; ++k )
	{
		//================ Pivotting ================//		
		amax = 0;
		mu = -1;
		lambda = -1;
		for( int i=k; i<m; ++i )
		{
			for( int j=k; j<n; ++j )
			{
				if( abs(a(i, j)) > amax )
				{
					amax = abs(a(i, j));
					mu = i;
					lambda = j;
				}
			}
		}

		// Swap row mu and k
		if( mu != -1 )
		{
			p.Swap( mu, k );
			a.SwapRow( mu, k );
		}
		
		// Swap column lambda and k
		if( lambda != -1 )
		{
			q.Swap( lambda, k );
			a.SwapColumn( lambda, k );
		}
		
		if( a(k, k) == 0 )
			continue;

		//============== LU decomposition ============//
		for( int i=k+1; i<m; ++i )
		{
			a(i, k) = a(i, k) / a(k, k);

			for( int j=k+1; j<n; ++j )
				a(i, j) -= a(i, k) * a(k, j);
		}

	}// end of k loop

	return true;
}



template< typename T >
inline bool LU_FullPivot( IMatrix<T>& lu, OreOreLib::Array<int>& p, OreOreLib::Array<int>& q, const IMatrix<T>& a )
{
	return lu.Copy(a) ? LU_FullPivot(lu, p, q) : false;
}





// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.330.2112&rep=rep1&type=pdf <-
// http://nkl.cc.u-tokyo.ac.jp/13n/SolverPrecond.pdf
// http://nkl.cc.u-tokyo.ac.jp/11e/05-Solver/solver.pdf


template< typename T >
inline bool ILU0( IMatrix<T>& a )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();

	for( int i=1; i<n; ++i )
	{
		for( int k=0; k<i; ++k )
		{
			if( a(i, k)==0 || a(k, k)==0 )
				continue;

			a(i, k) /= a(k, k);

			for( int j=k+1; j<n; ++j )
			{
				if( a(i, j) !=0 )
					a(i, j) -= a(i, k) * a(k, j);
			}
		}
	}

	return true;
}


template< typename T >
inline bool ILU0( IMatrix<T>&lu, const IMatrix<T>& a )
{
	return lu.Copy(a) ? ILU0(lu) : false;
}



// https://people.cs.kuleuven.be/~karl.meerbergen/didactiek/h03g1a/ilu.pdf
// https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs13160-017-0277-5/MediaObjects/13160_2017_277_Figa_HTML.gif <-
// Modified ILU(0)
template< typename T >
inline bool MILU0( IMatrix<T>& a )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();
	T s;

	for( int i=0; i<m; ++i )
	{
		s = 0;
		for( int k=0; k<i; ++k )
		{
			if( a(i, k)==0 || a(k, k)==0 )
				continue;

			a(i, k) /= a(k, k);

			for( int j=k+1; j<n; ++j )
			{
				if( a(k, j)==0 )
					continue;

				if( a(i, j)==0 )
					s += a(i, k) * a(k, j);
				else
					a(i, j) -= a(i, k) * a(k, j);

			}// end of j loop
		}// end of k loop

		a(i, i) -= s;

	}// end of i loop

	return true;
}


template< typename T >
inline bool MILU0( IMatrix<T>&lu, const IMatrix<T>& a )
{
	return lu.Copy(a) ? MILU0(lu) : false;
}




// https://www.jstage.jst.go.jp/article/journalam1998/7/0/7_0_279/_pdf
// https://fse.studenttheses.ub.rug.nl/11132/1/Koen_van_Geffen_2013_TWB.pdf
// https://www-users.cs.umn.edu/~saad/PDF/umsi-92-38.pdf

template< typename T >
inline bool ILUT( IMatrix<T>& l, IMatrix<T>& u, /*OreOreLib::Array<T>& w_, OreOreLib::Array<T>& d,*/ const IMatrix<T>& a, const T& t )
{
	if( !a.IsSquare() )
		return false;

	const int m = a.numRows();
	const int n = a.numCols();

	MatrixView<T> w(a, 0, 0, 1, n );

/*
	for( int i=1; i<n; ++i )
	{
//		for( int j=0; j<n; ++j )
//			w[j] = a(i, j);
		w.Init(a, i, 0, 1, n  );

		for( int k=0; k<i; ++k )
		{
			if( w(k)==0 ) continue;

			//tcout << "a(" << k << ", " << k << "): " << a(k, k) << tendl;
			//w(k) /= a(k, k);
			if( a(k,k)==0 ) w(k)=0;
			else w(k) /= a(k, k);
			w(k) *= T(w(k) >= t);// Apply drop rule to wk
			
			if( w(k) != 0 )
			{
				for( int j=k+1; j<n; ++j )
				{
					tcout << j << ":" << w(k)  <<tendl;
					w(j) -= w(k) * u(k, j);
					
				}
				//w.Display();
			}
		}
		tcout << tendl;
		// Apply drop rule to w
		for( int j=0; j<n; ++j )
			w(j) *= T(w(j) >= t);

		for( int j=0; j<i; ++j )
			l(i, j) = w(j);

		for( int j=i; j<n; ++j )
			u(i, j) = w(j);

		w.Clear();
		
//		l.Display();
//		u.Display();

	}

	return true;
*/	


	for( int i=0; i<n; ++i )
	{
		//T ti = 0;

		//for( int j=0; j<n; ++j )
		//{
			//w[j] = a(i, j);
			//ti += (T)pow( w[j], 2 );
		//}
		w.Init(a, i, 0, 1, n  );

		//ti = t * sqrt(ti);

		for( int k=0; k<i; ++k )
		{
			w(k) /= u(k, k);//d[k];
			
			if( abs(w(k)) >=t )
			{
				l(i, k) = w(k);
				//tcout << "l:" << i <<"," << k << tendl;
			}

			for( int j=k+1; j<n; ++j )
				if( u(k, j) != 0 )
					w(j) -= w(k) * u(k, j);
		}

		/*d[i]*/u(i, i) = w(i);

		for( int j=i+1; j<n; ++j )
			if( abs(w(j)) >= t )
				u(i, j) = w(j);

		l(i, i) = 1;
	}
	
	
	return true;
}



//template< typename T >
//inline bool ILUT( IMatrix<T>&l, IMatrix<T>&u, const IMatrix<T>& a, const T& t )
//{
//	int n = Min( a.numRows(), a.numCols() );
//	//OreOreLib::Array<T> w_(n);//, d(n);
//	return ILUT(l, u, /*w_, d,*/ a, t);
//}





//================== Inverse matrix calculation using LU decomposition ==============================//

// https://math.stackexchange.com/questions/1009916/easy-way-to-calculate-inverse-of-an-lu-decomposition/1013950
// A * A_inv = I
// L * U * A_inv = I
//
// for i=0,...,I.numCols()-1
//   Xi = X = U * A_inv(col[i])
//   Ii = I(col[i])
//   Solve L * Xi = Ii
//   Solve U * A_inv(col[i]) = Xi



// https://mobiusfunction.wordpress.com/2010/12/21/the-inverse-of-a-triangular-matrix-using-forward-substitution/
//  Step1. Calculate inverse matrices L^-1 with forwadd substitution, and U^-1 with backward substitution. 
//  Step2. A^-1 = U^-1 * L^-1
// Inverse_LU
template< typename T >
inline void Inverse_LU( IMatrix<T>& out, const IMatrix<T>& in )
{
	assert( in.numRows()==out.numRows() && in.numCols()==out.numCols() );

	const int m = in.numRows();
	const int n = in.numCols();

	DynamicMatrix<T> L_inv( m, n ), U_inv( m, n );

	// Use out as decomposed_lu temporarily.
	out = in;
	LU( out );

	// Calculate L_inv and U_inv
	for( int k=0; k<n; ++k )
	{
		// Solve L * L_inv = I with forward substitution
		for( int i=k; i<m; ++i )
		{
			L_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=0; j<i; ++j )
				L_inv(i, k) -= out(i, j) * L_inv(j, k);
		}

		// Solve U * U_inv = I with backward substitution
		for( int i=k; i>=0; i-- )
		{
			U_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=i+1; j<n; ++j )
				U_inv(i, k) -= out(i, j) * U_inv(j, k);
			U_inv(i, k) /= out(i, i);
		}
	}

	//L_inv.Display();
	//U_inv.Display();

	// in^-1 = U_inv * L_inv
	Multiply( out, U_inv, L_inv );
}



template< typename T >
inline void Inverse_LU_Pivot( IMatrix<T>& out, const IMatrix<T>& in )
{
	assert( in.numRows()==out.numRows() && in.numCols()==out.numCols() );

	const int m = in.numRows();
	const int n = in.numCols();

	DynamicMatrix<T> L_inv( m, n ), U_inv( m, n );
	OreOreLib::Array<int> pivot( m );

	// Use out as decomposed_lu temporarily.
	out = in;
	LU_Pivot( out, pivot );


	// Calculate L_inv and U_inv
	for( int k=0; k<n; ++k )
	{
		// Solve L * L_inv = I with forward substitution
		for( int i=k; i<m; ++i )
		{
			L_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=0; j<i; ++j )
				L_inv(i, k) -= out(i, j) * L_inv(j, k);
		}

		// Solve U * U_inv = I with backward substitution
		for( int i=k; i>=0; --i )
		{
			U_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=i+1; j<n; ++j )
				U_inv(i, k) -= out(i, j) * U_inv(j, k);
			U_inv(i, k) /= out(i, i);
		}
	}

	// in^-1 = U_inv * L_inv
	for( int i=0; i<m; ++i )
	{
		for( int j=0; j<n; ++j )
		{
			T& dest	= out(i, pivot[j]);// set outout element according to pivot
			dest	= 0;
			for( int k=0; k<n; ++k )
				dest += U_inv(i, k) * L_inv(k, j);
		}
	}

}



template< typename T >
inline void Inverse_LU_FullPivot( IMatrix<T>& out, const IMatrix<T>& in )
{
	assert( in.numRows()==out.numRows() && in.numCols()==out.numCols() );

	const int m = in.numRows();
	const int n = in.numCols();

	DynamicMatrix<T> L_inv( m, n ), U_inv( m, n );
	OreOreLib::Array<int> pivot( m ), colpivot( n );

	// Use out as decomposed_lu temporarily.
	out = in;
	LU_FullPivot( out, pivot, colpivot );


	// Calculate L_inv and U_inv
	for( int k=0; k<n; ++k )
	{
		// Solve L * L_inv = I with forward substitution
		for( int i=k; i<m; ++i )
		{
			L_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=0; j<i; ++j )
				L_inv(i, k) -= out(i, j) * L_inv(j, k);
		}

		// Solve U * U_inv = I with backward substitution
		for( int i=k; i>=0; --i )
		{
			U_inv(i, k) = T(i==k);// (i, k) element of idendity matrix I
			for( int j=i+1; j<n; ++j )
				U_inv(i, k) -= out(i, j) * U_inv(j, k);
			U_inv(i, k) /= out(i, i);
		}
	}

	// in^-1 = U_inv * L_inv
	for( int i=0; i<m; ++i )
	{
		for( int j=0; j<n; ++j )
		{
			T& dest	= out(colpivot[i], pivot[j]);// set outout element according to pivot
			dest	= 0;
			for( int k=0; k<n; ++k )
				dest += U_inv(i, k) * L_inv(k, j);
		}
	}

}




template< typename T >
class LU_Solver : public Solver<T>
{
public:

	LU_Solver()
	{
	}

	LU_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~LU_Solver()
	{
		m_Comp.Release();
		m_y.Release();

		m_Inverse.Release();
		L_inv.Release();
		U_inv.Release();
		m_bInvUpdated = false;
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_Comp.Init( m, n );
		m_y.Init( n, batchSize );

		L_inv.Init( m, n );
		U_inv.Init( m, n );
		m_Inverse.Init( m, n );

		SetCoeff( A );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		m_bInvUpdated = false;
		return LU( m_Comp, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Comp.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Comp.numCols(), batchSize );
			return true;
		}
		return false;
	}


	//void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_Comp.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(i, 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_Comp(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
	//	for( int i=m_Comp.numRows()-1; i>=0; --i )
	//	{
	//		x(i, 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_Comp.numCols(); ++j )
	//			x(i, 0) -= m_Comp(i, j) * x(j, 0);
	//		x(i, 0) /= m_Comp(i, i);
	//	}
	//}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();
			
			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
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
			int tid = omp_get_thread_num();

			//tcout << k << tendl;
			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
			}
		}
	}

	#endif


	const DynamicMatrix<T>& Inverse()
	{
		if(m_bInvUpdated==true)
			return m_Inverse;

		// Calculate L_inv and U_inv
		for( int n=0; n<L_inv.numCols(); ++n )
		{
			// Solve L * L_inv = I with forward substitution
			for( int i=n; i<m_Comp.numRows(); ++i )
			{
				L_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=0; j<i; ++j )
					L_inv(i, n) -= m_Comp(i, j) * L_inv(j, n);
			}

			// Solve U * U_inv = I with backward substitution
			for( int i=n; i>=0; --i )
			{
				U_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					U_inv(i, n) -= m_Comp(i, j) * U_inv(j, n);
				U_inv(i, n) /= m_Comp(i, i);
			}
		}

		// in^-1 = U_inv * L_inv
		for( int i=0; i<m_Inverse.numRows(); ++i )
		{
			for( int j=0; j<m_Inverse.numCols(); ++j )
			{
				T& dest	= m_Inverse(i, j);
				dest	= 0;
				for( int k=0; k<U_inv.numCols(); ++k )
					dest += U_inv(i, k) * L_inv(k, j);
			}
		}

		m_bInvUpdated = true;

		return m_Inverse;
	}


	const T Determinant()
	{
		T det = 0;
		for( int i=0; i<m_Comp.numRows(); ++i )
			det += pow( m_Comp(i, i), 2 );
		return det;
	}



private:

	DynamicMatrix<T>				m_Comp;	// Composed LU matrix
	DynamicMatrix<T>				m_y;	// Intermediate variables.

	// temporal buffers for invers matrix calculation
	DynamicMatrix<T>				m_Inverse;
	DynamicMatrix<T>				L_inv, U_inv;
	bool					m_bInvUpdated;
};




template< typename T >
class PivotLU_Solver : public Solver<T>
{
public:

	PivotLU_Solver()
	{
	}


	PivotLU_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~PivotLU_Solver()
	{
		m_Comp.Release();
		m_y.Release();
		m_Pivot.Release();

		m_Inverse.Release();
		L_inv.Release();
		U_inv.Release();
		m_bInvUpdated = false;
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_Comp.Init( m, n );
		m_y.Init( n, batchSize );
		m_Pivot.Init( n );

		L_inv.Init( m, n );
		U_inv.Init( m, n );
		m_Inverse.Init( m, n );

		SetCoeff( A );

		//for( int i=0; i<m_Pivot.Length(); ++i )
		//	tcout << m_Pivot[i] << tendl;
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		m_bInvUpdated = false;
		return LU_Pivot( m_Comp, m_Pivot, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Comp.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Comp.numCols(), batchSize );
			return true;
		}
		return false;
	}


	//void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_Comp.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(m_Pivot[i], 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_Comp(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
	//	for( int i=m_Comp.numRows()-1; i>=0; --i )
	//	{
	//		if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.
	//
	//		x(i, 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_Comp.numCols(); ++j )
	//			x(i, 0) -= m_Comp(i, j) * x(j, 0);
	//		x(i, 0) /= m_Comp(i, i);
	//	}
	//}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(m_Pivot[i], k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
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
			int tid = omp_get_thread_num();

			//tcout << k << tendl;
			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(m_Pivot[i], k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
			}
		}
	}

	#endif


	const DynamicMatrix<T>& Inverse()
	{
		if(m_bInvUpdated==true)
			return m_Inverse;

		// Calculate L_inv and U_inv
		for( int n=0; n<L_inv.numCols(); ++n )
		{
			// Solve L * L_inv = I with forward substitution
			for( int i=n; i<m_Comp.numRows(); ++i )
			{
				L_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=0; j<i; ++j )
					L_inv(i, n) -= m_Comp(i, j) * L_inv(j, n);
			}

			// Solve U * U_inv = I with backward substitution
			for( int i=n; i>=0; --i )
			{
				U_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					U_inv(i, n) -= m_Comp(i, j) * U_inv(j, n);
				U_inv(i, n) /= m_Comp(i, i);
			}
		}

		// in^-1 = U_inv * L_inv
		for( int i=0; i<m_Inverse.numRows(); ++i )
		{
			for( int j=0; j<m_Inverse.numCols(); ++j )
			{
				T& dest	= m_Inverse(i, m_Pivot[j]);// set outout element according to pivot
				dest	= 0;
				for( int k=0; k<U_inv.numCols(); ++k )
					dest += U_inv(i, k) * L_inv(k, j);
			}
		}

		m_bInvUpdated = true;

		return m_Inverse;
	}


	const T Determinant()
	{
		T det = 0;
		for( int i=0; i<m_Comp.numRows(); ++i )
			det += pow( m_Comp(i, i), 2 );
		return det;
	}



private:

	DynamicMatrix<T>				m_Comp;	// Composed LU matrix
	DynamicMatrix<T>				m_y;	// Intermediate variables.
	OreOreLib::Array<int>	m_Pivot;// Variable order for pivoting
	
	// temporal buffers for invers matrix calculation
	DynamicMatrix<T>				m_Inverse;
	DynamicMatrix<T>				L_inv, U_inv;
	bool					m_bInvUpdated;
};




template< typename T >
class FullPivotLU_Solver : public Solver<T>
{
public:

	FullPivotLU_Solver()
	{
	}


	FullPivotLU_Solver( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~FullPivotLU_Solver()
	{
		m_Comp.Release();
		m_y.Release();
		m_RowPivot.Release();
		m_ColPivot.Release();

		m_Inverse.Release();
		L_inv.Release();
		U_inv.Release();
		m_bInvUpdated = false;
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		int m = A.numRows();
		int n = A.numCols();

		m_Comp.Init( m, n );
		m_y.Init( n, batchSize );
		m_RowPivot.Init( m );
		m_ColPivot.Init( n );

		L_inv.Init( m, n );
		U_inv.Init( m, n );
		m_Inverse.Init( m, n );

		SetCoeff( A );

		//for( int i=0; i<m_Pivot.Length(); ++i )
		//	tcout << m_Pivot[i] << tendl;
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		m_bInvUpdated = false;
		return LU_FullPivot( m_Comp, m_RowPivot, m_ColPivot, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Comp.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Comp.numCols(), batchSize );
			return true;
		}
		return false;
	}


	//void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_Comp.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(m_Pivot[i], 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_Comp(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
	//	for( int i=m_Comp.numRows()-1; i>=0; --i )
	//	{
	//		if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.
	//
	//		x((m_ColPivot[i], 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_Comp.numCols(); ++j )
	//			x((m_ColPivot[i], 0) -= m_Comp(i, j) * x(m_ColPivot[j], 0);
	//		x((m_ColPivot[i], 0) /= m_Comp(i, i);
	//	}
	//}


	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(m_RowPivot[i], k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(m_ColPivot[i], k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(m_ColPivot[i], k) -= m_Comp(i, j) * x(m_ColPivot[j], k);
				x(m_ColPivot[i], k) /= m_Comp(i, i);
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
			int tid = omp_get_thread_num();

			//tcout << k << tendl;
			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(m_RowPivot[i], k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(m_ColPivot[i], k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(m_ColPivot[i], k) -= m_Comp(i, j) * x(m_ColPivot[j], k);
				x(m_ColPivot[i], k) /= m_Comp(i, i);
			}
		}
	}

	#endif


	const DynamicMatrix<T>& Inverse()
	{
		if(m_bInvUpdated==true)
			return m_Inverse;

		// Calculate L_inv and U_inv
		for( int n=0; n<L_inv.numCols(); ++n )
		{
			// Solve L * L_inv = I with forward substitution
			for( int i=n; i<m_Comp.numRows(); ++i )
			{
				L_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=0; j<i; ++j )
					L_inv(i, n) -= m_Comp(i, j) * L_inv(j, n);
			}

			// Solve U * U_inv = I with backward substitution
			for( int i=n; i>=0; --i )
			{
				U_inv(i, n) = T(i==n);// (i, n) element of idendity matrix I
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					U_inv(i, n) -= m_Comp(i, j) * U_inv(j, n);
				U_inv(i, n) /= m_Comp(i, i);
			}
		}

		// in^-1 = U_inv * L_inv
		for( int i=0; i<m_Inverse.numRows(); ++i )
		{
			for( int j=0; j<m_Inverse.numCols(); ++j )
			{
				T& dest	= m_Inverse(m_ColPivot[i], m_RowPivot[j]);// set outout element according to pivot
				dest	= 0;
				for( int k=0; k<U_inv.numCols(); ++k )
					dest += U_inv(i, k) * L_inv(k, j);
			}
		}

		m_bInvUpdated = true;

		return m_Inverse;
	}


	const T Determinant()
	{
		T det = 0;
		for( int i=0; i<m_Comp.numRows(); ++i )
			det += pow( m_Comp(i, i), 2 );
		return det;
	}



private:

	DynamicMatrix<T>				m_Comp;	// Composed LU matrix
	DynamicMatrix<T>				m_y;	// Intermediate variables.
	OreOreLib::Array<int>	m_RowPivot, m_ColPivot;// Variable order for pivoting
	
	// temporal buffers for invers matrix calculation
	DynamicMatrix<T>				m_Inverse;
	DynamicMatrix<T>				L_inv, U_inv;
	bool					m_bInvUpdated;
};




template< typename T >
class ILU0_Preconditioner : public Solver<T>
{
public:

	ILU0_Preconditioner()
	{
	}


	ILU0_Preconditioner( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~ILU0_Preconditioner()
	{
		m_Comp.Release();
		m_y.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		m_Comp.Init( A.numRows(), A.numCols() );
		m_y.Init( A.numCols(), batchSize );

		SetCoeff( A );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		return ILU0( m_Comp, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Comp.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Comp.numCols(), batchSize );
			return true;
		}
		return false;
	}


	//void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_Comp.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(i, 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_Comp(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
	//	for( int i=m_Comp.numRows()-1; i>=0; --i )
	//	{
	//		if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.
	//
	//		x(i, 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_Comp.numCols(); ++j )
	//			x(i, 0) -= m_Comp(i, j) * x(j, 0);
	//		x(i, 0) /= m_Comp(i, i);
	//	}
	//}

	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
			}
		}
	}



private:

	DynamicMatrix<T>	m_Comp;	// Composed LU matrix
	DynamicMatrix<T>	m_y;	// Intermediate variables.

};




// Modified ILU(0) Preconditioner
template< typename T >
class MILU0_Preconditioner : public Solver<T>
{
public:

	MILU0_Preconditioner()
	{
	}


	MILU0_Preconditioner( const IMatrix<T>& A, int batchSize=1 )
	{
		Init( A, batchSize );
	}


	~MILU0_Preconditioner()
	{
		m_Comp.Release();
		m_y.Release();
	}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 )
	{
		m_Comp.Init( A.numRows(), A.numCols() );
		m_y.Init( A.numCols(), batchSize );

		SetCoeff( A );
	}


	virtual bool SetCoeff( const IMatrix<T>& A )
	{
		return MILU0( m_Comp, A );
	}


	virtual bool SetBatchSize( int batchSize )
	{
		assert( batchSize>0 );
		if( !m_Comp.IsEmpty() && m_y.numCols()!=batchSize )
		{
			m_y.Init( m_Comp.numCols(), batchSize );
			return true;
		}
		return false;
	}


	//virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_Comp.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(i, 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_Comp(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
	//	for( int i=m_Comp.numRows()-1; i>=0; --i )
	//	{
	//		if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.
	//
	//		x(i, 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_Comp.numCols(); ++j )
	//			x(i, 0) -= m_Comp(i, j) * x(j, 0);
	//		x(i, 0) /= m_Comp(i, i);
	//	}
	//}

	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_Comp.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Forward substitution using lower-left part
			for( int i=0; i<m_Comp.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_Comp(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part// BackwardSubstitution( m_Comp, x, m_y );
			for( int i=m_Comp.numRows()-1; i>=0; --i )
			{
				if( m_Comp(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_Comp.numCols(); ++j )
					x(i, k) -= m_Comp(i, j) * x(j, k);
				x(i, k) /= m_Comp(i, i);
			}
		}
	}



private:

	DynamicMatrix<T>	m_Comp;	// Composed LU matrix
	DynamicMatrix<T>	m_y;	// Intermediate variables.

};




template< typename T >
class ILUT_Preconditioner : public Solver<T>
{
public:

	ILUT_Preconditioner()
	{
	}


	ILUT_Preconditioner( const IMatrix<T>& A, const T& t, int batchSize=1 )
	{
		Init( A, t, batchSize );
	}


	~ILUT_Preconditioner()
	{
		m_L.Release();
		m_U.Release();
		m_y.Release();
	}


	virtual void Init( const IMatrix<T>& A, const T& t, int batchSize=1 )
	{
		m_L.Init( A.numRows(), A.numCols() );
		m_U.Init( A.numRows(), A.numCols() );
		m_y.Init( A.numCols(), batchSize );

		SetCoeff( A, t );
	}


	bool SetCoeff( const IMatrix<T>& A, const T& t )
	{
		return ILUT( m_L, m_U, A, t );
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


	//virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	//{
	//	if( !(x.IsSameSize(b) && x.numRows()==m_L.numRows()) )
	//		return;

	//	// Forward substitution using lower-left part
	//	for( int i=0; i<m_L.numRows(); ++i )
	//	{
	//		m_y(i, 0) = b(i, 0);
	//		for( int j=0; j<i; ++j )
	//			m_y(i, 0) -= m_L(i, j) * m_y(j, 0);
	//	}

	//	// Backward substitution using upper-right part// BackwardSubstitution( m_U, x, m_y );
	//	for( int i=m_U.numRows()-1; i>=0; --i )
	//	{
	//		if( m_U(i, i)==0 )	continue;// Ignore ZERO DIVISION.
	//
	//		x(i, 0) = m_y(i, 0);
	//		for( int j=i+1; j<m_U.numCols(); ++j )
	//			x(i, 0) -= m_U(i, j) * x(j, 0);
	//		x(i, 0) /= m_U(i, i);
	//	}
	//}

	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )
	{
		if( !(x.IsSameSize(b) && x.numRows()==m_L.numRows()) )
			return;

		for( int k=0; k<b.numCols(); ++k )
		{
			int tid = k % m_y.numCols();

			// Forward substitution using lower-left part
			for( int i=0; i<m_L.numRows(); ++i )
			{
				m_y(i, tid) = b(i, k);
				for( int j=0; j<i; ++j )
					m_y(i, tid) -= m_L(i, j) * m_y(j, tid);
			}

			// Backward substitution using upper-right part// BackwardSubstitution( m_U, x, m_y );
			for( int i=m_U.numRows()-1; i>=0; --i )
			{
				if( m_U(i, i)==0 )	continue;// Ignore ZERO DIVISION.

				x(i, k) = m_y(i, tid);
				for( int j=i+1; j<m_U.numCols(); ++j )
					x(i, k) -= m_U(i, j) * x(j, k);
				x(i, k) /= m_U(i, i);
			}
		}
	}



private:

	DynamicMatrix<T>	m_L, m_U;
	DynamicMatrix<T>	m_y;	// Intermediate variables.


	bool SetCoeff( const IMatrix<T>& A ){ return false; }


};






//##############################################################################//
//							Experimental implementation							//
//##############################################################################//


//template< typename T >
//void LUCrout( IMatrix<T>& L, IMatrix<T>& U, const IMatrix<T>& A )
//{
//	for( int i=0; i<A.numRows(); ++i )
//	{
//		U(i, i) = 1;
//		// l(i,j)
//		for( int j=0; j<=i; ++j )
//		{
//			T l_ij = A(i, j);
//			for( int k=0; k<j; ++k )
//				l_ij -= L(i, k) * U(k, j);
//
//			L(i, j) = l_ij;
//		}
//		
//		// u(i,j)
//		for( int j=i+1; j<A.numCols(); ++j )
//		{
//			T u_ij = A(i, j);
//			for( int k=0; k<i; ++k )
//				u_ij -= L(i, k) * U(k, j);
//
//			U(i, j) = u_ij / L(i, i);
//		}
//
//	}// end of i loop
//
//}




//template< typename T >
//void LUDoolittle( IMatrix<T>& L, IMatrix<T>& U, const IMatrix<T>& A )
//{
//	for( int i=0; i<A.numRows(); ++i )
//	{
//		L(i, i) = 1;
//		// l(i,j)
//		for( int j=0; j<i; ++j )
//		{
//			T l_ij = A(i, j);
//			for( int k=0; k<j; ++k )
//				l_ij -= L(i, k) * U(k, j);
//
//			L(i, j) = l_ij / U(j, j);
//		}
//
//		// u(i,j)
//		for( int j=i; j<A.numCols(); ++j )
//		{
//			T u_ij = A(i, j);
//			for( int k=0; k<i; ++k )
//				u_ij -= L(i, k) * U(k, j);
//
//			U(i, j) = u_ij;
//		}
//
//	}// end of i loop
//
//}



#endif // !LU_H
