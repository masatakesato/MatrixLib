#ifndef I_MATRIX_H
#define	I_MATRIX_H

#include	<oreore/common/Utility.h>
#include	<oreore/common/TString.h>
#include	<oreore/mathlib/MathLib.h>



// memory acces method
enum MEM
{
	PARTIAL,// Using partial area of allocated memries.
	ENTIRE	// Using entire allocated memory.
};



const int DynamicDim = -1;



template< typename T >
class IMatrix
{
public:

	// Default constructor
	IMatrix( MEM mem=MEM::PARTIAL )
		: m_numRows( 0 )
		, m_numCols( 0 )
		, memcopiable( mem )
	{
	}


	// Constructor
	IMatrix( int numrows, int numcols, MEM mem )
		: m_numRows( numrows )
		, m_numCols( numcols )
		, memcopiable( mem )
	{
		assert( numrows>0 && numcols>0 );
	}


	// Copy constructor
	//IMatrix( const IMatrix& obj )
	//	: m_numRows( obj.m_numRows )
	//	, m_numCols( obj.m_numCols )
	//	, memcopiable( obj.memcopiable )
	//{
	//	
	//}


	// Move constructor
	//IMatrix( IMatrix&& obj )
	//	: m_numRows( obj.m_numRows )
	//	, m_numCols( obj.m_numCols )
	//	, memcopiable( obj.memcopiable )
	//{
	//	
	//}


	virtual ~IMatrix()
	{

	}



	// Copy assignment operator. Derived class MUST have implementation.
	virtual IMatrix& operator=( const IMatrix& obj )
	{
		m_numRows	= obj.m_numRows;
		m_numCols	= obj.m_numCols;
		
		//memcopiable = obj.memcopiable;// Derived class must handle this precedure
		//Copy( obj );// Derived class must handle this precedure

		return *this;
	}


	// Move assignment operator. Derived class MUST have implementation.
	virtual IMatrix& operator=( IMatrix&& obj )
	{
		m_numRows = obj.m_numRows;
		m_numCols = obj.m_numCols;

		obj.m_numRows = 0;
		obj.m_numCols = 0;
		//memcopiable = obj.memcopiable;// Derived class must handle this precedure

		return *this;
	}



	// Compound assignment operator +=
	template< typename MatType >
	IMatrix& operator+=( const IMatrix<MatType>& right )
	{
		assert( IsSameSize(right) );

		for( int i=0; i<m_numRows * m_numCols; ++i )
			(*this)(i) += right(i);

		return *this;
	}


	// Compound assignment operator -=
	template< typename MatType >
	IMatrix& operator-=( const IMatrix<MatType>& right )
	{
		assert( IsSameSize(right) );

		for( int i=0; i<m_numRows * m_numCols; ++i )
			(*this)(i) -= right(i);

		return *this;
	}


	//// Compound assignment operator *=. Derived class MUST have implementation.
	//template< typename MatType >
	//IMatrix& operator*=( const IMatrix<MatType>& right )
	//{
	//	return static_cast<Derived*>(this)->operator*=( static_cast<const Derived&>(right) );
	//}


	virtual void Clear()
	{
		memset( pData(), 0, Size() );
	}


	virtual bool Copy( const IMatrix& src )
	{
		if( !src.IsSameSize(*this) )
			return false;

		for( int i=0; i<numElms(); ++i )
			(*this)(i) = src(i);

		return true;
	}


	virtual bool CopyTo( IMatrix& dst ) const
	{
		if( !dst.IsSameSize(*this) )
			return false;

		for( int i=0; i<dst.numElms(); ++i )
			dst(i) = (*this)(i);

		return true;
	}


	virtual bool CopyRow( const IMatrix& src, int row_src, int row_target )
	{
		if( m_numCols!=src.m_numCols || src.m_numRows<=row_src || m_numRows<=row_target )
			return false;

		//for( int i=0; i<m_numCols; ++i )
		//	(*this)(row_target, i) -= src(row_src, i);
		memcpy_s( ptr(row_target, 0), m_numCols * sizeof(T), src.ptr(row_src, 0), src.m_numCols * sizeof(T) );

		return true;
	}


	virtual bool CopyRow( const T src[], int row_target )
	{
		if( m_numRows<=row_target )
			return false;

		//for( int i=0; i<m_numCols; ++i )
		//	(*this)(row_target, i) -= src[i];
		memcpy_s( ptr(row_target, 0), m_numCols * sizeof(T), src, m_numCols * sizeof(T) );

		return true;
	}


	virtual bool CopyColumn( const IMatrix& src, int col_src, int col_target )
	{
		if( m_numRows!=src.m_numRows || src.m_numCols<col_src || m_numCols<=col_target )
			return false;

		for( int i=0; i<m_numRows; ++i )
			(*this)(i, col_target) = src(i, col_src);

		return true;
	}


	virtual bool CopyColumn( const T src[], int col_target )
	{
		if( m_numCols<=col_target )
			return false;

		for( int i=0; i<m_numRows; ++i )
			(*this)(i, col_target) = src[i];

		return true;
	}


	void SwapRow( int i, int j )
	{
		assert( i>=0 && j>=0 && i<m_numRows && j<m_numRows );

		if( i==j ) return;

		for( int k=0; k<m_numCols; ++k )
		{
			T temp = (*this)(i, k);
			(*this)(i, k) = (*this)(j, k);
			(*this)(j, k) = temp;
		}
	}


	void SwapColumn( int i, int j )
	{
		assert( i>=0 && j>=0 && i<m_numCols && j<m_numCols );

		if( i==j ) return;

		for( int k=0; k<m_numRows; ++k )
		{
			T temp = (*this)(k, i);
			(*this)(k, i) = (*this)(k, j);
			(*this)(k, j) = temp;
		}
	}


	virtual void Transpose()
	{
		// Derived class should implement
	}


	bool IsFinite() const
	{
		for( int i=0; i<m_numRows*m_numCols; ++i )
			if( !isfinite( (*this)(i) ) ) return false;
		return !IsEmpty();
	}


	template< typename U >
	bool IsSameSize( const IMatrix<U>& m ) const
	{
		return m_numRows==/*m.numRows()*/m.m_numRows && m_numCols==/*m.numCols()*/m.m_numCols;
	}


	bool IsEmpty() const
	{
		return m_numRows*m_numCols==0;
	}


	bool IsSquare() const
	{
		return (m_numRows==m_numCols) && (m_numRows>0);
	}


	bool IsSymmetric() const
	{
		if( m_numRows != m_numCols )
			return false;

		for( int i=1; i<m_numRows-1; ++i )
			for( int j=i+1; j<m_numCols; ++j )
				if( (*this)(i, j) != (*this)(j, i) ) return false;

		return true;
	}


	int Size() const
	{
		return m_numRows * m_numCols * sizeof(T);
	}


	int numElms() const
	{
		return m_numRows * m_numCols;
	}


	int numRows() const
	{
		return m_numRows;
	}


	int numCols() const
	{
		return m_numCols;
	}


	void Display() const
	{
		tcout << typeid(*this).name() << _T( "[" << m_numRows << "][" << m_numCols << "]:\n" );
		for(int row=0; row<m_numRows; ++row)
		{
			tcout << " ";
			for(int col=0; col<m_numCols; ++col)
			{
				tcout << (*this)( row, col ) << _T( " " );
			}// end of col loop
			tcout << tendl;
		}// end of row loop
	}



	// Subscription operator for read only.( called if Matrix is const )
	const T& operator()( int i ) const&				{ return *ptr(i); }
	const T& operator()( int row, int col ) const&	{ return *ptr(row, col); }

	// Subscription operator for read-write.( called if Matrix is non-const )
	T& operator()( int i ) &						{ return *ptr(i); }
	T& operator()( int row, int col ) &				{ return *ptr(row, col); }

	// Subscription operator. ( called by following cases: "T& a = Matrix<T>(10,10)(n)", "auto&& a = Matrix<T>(20,4)(n)" )
	T operator()( int i ) const&&					{ return std::move(*ptr(i)); }
	T operator()( int row, int col ) const&&		{ return std::move(*ptr(row, col)); }


	//T *m( int x, int y ) const						{ return ptr(y, x); }// x: horizontal, y: vertical
	//T &at( int row, int col ) const					{ return *ptr(row, col); }

	template< typename Type=T >
	Type *pDataAs() const 							{ return (Type *)pData(); }

	template< typename Type=T >
	Type &DataAs( int i, int j ) const 				{ return (Type *)ptr(i, j); }


	//=========================== virtual methods =============================//
	virtual T *ptr( int row, int col ) const=0;//{ return nullptr; }// // row: vertical, col: horizontal
	virtual T *ptr( int i ) const=0;//{ return nullptr; }//
	virtual unsigned char *pData() const=0;//{ return nullptr; }//



protected:

	int	m_numRows;
	int	m_numCols;



public:

	const MEM memcopiable;

};



template< typename T, int R, int C, typename enable=void > class Matrix; 



#endif // !I_MATRIX_H
