#ifndef DYNAMIC_MATRIX_H
#define	DYNAMIC_MATRIX_H

#include	<oreore/memory/Memory.h>

#include	"IMatrix.h"




//##########################################################################//
//							N x M Matrix Class								//
//##########################################################################//


template< typename T >
using DynamicMatrix = Matrix< T, DynamicDim, DynamicDim >;




template< typename T >
using RowVector = Matrix< T, 1, DynamicDim >;


template< typename T >
using ColVector = Matrix< T, DynamicDim, 1 >;




template< typename T, int R, int C >
class Matrix< T, R, C, std::enable_if_t< std::is_arithmetic_v<T> && (R==DynamicDim || C==DynamicDim) > > : public IMatrix<T>
{
public:

	// Default constructor
	Matrix()
		: IMatrix<T>( MEM::ENTIRE )
		, m_pData( nullptr )
	{
	}


	// Constructor for full dynamic matrix
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	Matrix( int numrows, int numcols, const T* const pdata=nullptr )
		: IMatrix<T>(numrows, numcols, MEM::ENTIRE )
		, m_pData( new T[ this->numElms() ]() )
	{
		if( pdata )
			memcpy_s( m_pData, this->Size(), pdata, this->Size() );
	}


	// Constructor for row vector
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_==DynamicDim>* = nullptr >
	Matrix( int numcols, const T* const pdata=nullptr )
		: IMatrix<T>( R_, numcols, MEM::ENTIRE )
		, m_pData( new T[ this->numElms() ]() )
	{
		if( pdata )
			memcpy_s( m_pData, this->Size(), pdata, this->Size() );
	}


	// Constructor for col vector
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==1>* = nullptr >
	Matrix( int numrows, const T* const pdata=nullptr )
		: IMatrix<T>( numrows, C_, MEM::ENTIRE )
		, m_pData( new T[ this->numElms() ]() )
	{
		if( pdata )
			memcpy_s( m_pData, this->Size(), pdata, this->Size() );
	}


	// Constructor
	Matrix( const IMatrix<T>& obj )
		: IMatrix<T>( obj.numRows(), obj.numCols(), MEM::ENTIRE )
		, m_pData( nullptr )
	{
		if( obj.pData() )
		{
			m_pData = new T[ this->numElms() ];
			Copy( obj );
		}
	}


	// Constructor
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	Matrix( std::initializer_list<std::initializer_list<T>> ilist )
		: IMatrix<T>( MEM::ENTIRE )
		, m_pData( nullptr )
	{
		// Extract numRows and numCols from ilist
		this->m_numRows = int( ilist.size() );
		for( auto &row_data : ilist )
			this->m_numCols = Max( this->m_numCols, (int)row_data.size() );

		// Allocate memory
		m_pData	= new T[ this->numElms() ];

		// copy value
		auto p = m_pData;
		for( auto& row_data : ilist )
			for( const auto& val : row_data )
				*p++ = val;
	}


	// Constructor for row vector
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_==DynamicDim>* = nullptr >
	Matrix( std::initializer_list<T> ilist )
		: IMatrix<T>( R_, int(ilist.size()), MEM::ENTIRE )
		, m_pData( new T[ this->numElms() ]() )
	{
		auto p = m_pData;
		for( const auto& val : ilist )	*p++ = val;
	}


	// Constructor for col vector
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==1>* = nullptr >
	Matrix( std::initializer_list<T> ilist )
		: IMatrix<T>( int(ilist.size()), C_, MEM::ENTIRE )
		, m_pData( new T[ this->numElms() ]() )
	{
		auto p = m_pData;
		for( const auto& val : ilist )	*p++ = val;
	}


	// Copy constructor
	Matrix( const Matrix& obj )
		: IMatrix<T>( obj.m_numRows, obj.m_numCols, MEM::ENTIRE )
		, m_pData( nullptr )
	{
		if( obj.m_pData )
		{
			m_pData = new T[ this->numElms() ];
			memcpy_s( m_pData, this->Size(), obj.m_pData, obj.Size() );
		}
	}


	// Move constructor
	Matrix( Matrix&& obj )
		: IMatrix<T>( obj.m_numRows, obj.m_numCols, MEM::ENTIRE )
		, m_pData( obj.m_pData )
	{
		obj.m_numRows	= 0;
		obj.m_numCols	= 0;
		obj.m_pData		= nullptr;// clear reference from obj
	}


	virtual ~Matrix()
	{
		SafeDeleteArray( m_pData );
	}


	// Copy assignment operator
	Matrix& operator=( const Matrix& obj )
	{
		if( this != &obj )
		{
			if( !obj.IsSameSize(*this) )
				Allocate( obj.m_numRows, obj.m_numCols );

			memcpy_s( m_pData, this->Size(), obj.m_pData, obj.Size() );
		}

		return *this;
	}



	Matrix& operator=( const IMatrix<T>& obj )
	{
		if( this != &obj )
		{
			if( !obj.IsSameSize(*this) )
				Allocate( obj.numRows(), obj.numCols() );

			Copy(obj);
		}

		return *this;
	}


	// Move assignment operator
	Matrix& operator=( Matrix&& obj )
	{
		if( this != &obj )
		{
			// free current m_pData first.
			SafeDeleteArray( m_pData );

			// copy data to *this
			this->m_numRows	= obj.m_numRows;
			this->m_numCols	= obj.m_numCols;

			m_pData		= obj.m_pData;

			// clear reference from src
			obj.m_numRows	= 0;
			obj.m_numCols	= 0;
			obj.m_pData	= nullptr;
		}

		return *this;
	}


	// Matrix-specific implementation
	template< typename U >
	std::enable_if_t<std::is_arithmetic_v<U>, Matrix& >
	operator*=( const IMatrix<U>& right )
	{
		assert( this->m_numCols==right.numRows() );

		DynamicMatrix<U> buf( this->m_numRows, right.numCols() );

		for( int i=0; i<buf.m_numRows; ++i )
		{
			for( int j=0; j<buf.m_numCols; ++j )
			{
				for( int k=0; k<this->m_numCols; ++k )
					buf(i, j) += (*this)(i, k) * right(k, j);
			}
		}

		if( !buf.IsSameSize(*this) )
			Allocate( buf.m_numRows, buf.m_numCols );
		Copy( buf );

		//*this = std::move(buf);// TODO: MatrixViewの参照先が無効化する. ハンドル実装してからムーブ使う

		return *this;
	}


	// Init full dynamic matrix
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==DynamicDim && C_==DynamicDim, void >
	Init( int numrows, int numcols )
	{
		Allocate( numrows, numcols );
	}


	// Init row vector
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==1 && C_==DynamicDim, void >
	Init( int numcols )
	{
		Allocate( R_, numcols );
	}


	// Init col vector
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==DynamicDim && C_==1, void >
	Init( int numrows )
	{
		Allocate( numrows, C_ );
	}


	void Release()
	{
		this->m_numRows	= 0;
		this->m_numCols	= 0;

		SafeDeleteArray( m_pData );
	}


	void Clear()
	{
		memset( m_pData, 0, this->Size() );
	}


	void Resize( int numrows, int numcols )
	{
		assert( numrows>0 && numcols>0 );

		Release();

		this->m_numRows	= numrows;
		this->m_numCols	= numcols;

		m_pData			= new T[ numrows * numcols ]();

	}



	bool Copy( const IMatrix<T>& src )
	{
		if( !src.IsSameSize(*this) )
			return false;

		if( !src.memcopiable )
			return src.CopyTo( *this );

		memcpy_s( pData(), this->Size(), src.pData(), src.Size() );
		return true;
	}


	bool CopyTo( IMatrix<T>& dst ) const
	{
		if( !dst.IsSameSize(*this) )
			return false;

		if( !dst.memcopiable )
			return dst.Copy( *this );

		memcpy_s( dst.pData(), dst.Size(), pData(), this->Size() );
		return true;
	}


	// Same as IMatrix implementation
	//bool CopyRow( const Matrix<T>& src, int row_src, int row_target )
	//{
	//	if( this->m_numCols!=src.m_numCols || src.m_numRows<=row_src || this->m_numRows<=row_target )
	//		return false;

	//	memcpy_s( ptr(row_target, 0), this->m_numCols * sizeof(T), src.ptr(row_src, 0), src.m_numCols * sizeof(T) );
	//	return true;
	//}


	// Same as IMatrix implementation
	//bool CopyRow( const T src[], int row_target )
	//{
	//	if( this->m_numRows<=row_target )
	//		return false;

	//	memcpy_s( ptr(row_target, 0), this->m_numCols * sizeof(T), src, this->m_numCols * sizeof(T) );
	//	return true;
	//}


	// Same as IMatrix implementation
	//bool CopyColumn( const Matrix<T>& src, int col_src, int col_target )
	//{
	//	if( this->m_numRows!=src.m_numRows || src.m_numCols<col_src || this->m_numCols<=col_target )
	//		return false;

	//	for( int i=0; i<this->m_numRows; ++i )
	//		(*this)(i, col_target) = src(i, col_src);

	//	return true;
	//}


	// Same as IMatrix implementation
	//bool CopyColumn( const T src[], int col_target )
	//{
	//	if( m_numCols<=col_target )
	//		return false;

	//	for( int i=0; i<m_numRows; ++i )
	//		(*this)(i, col_target) = src[i];

	//	return true;
	//}


	void Transpose()
	{
		T *temp	= new T[ this->numElms() ];
		memcpy_s( temp, this->Size(), pData(), this->Size() );

		int r = this->m_numRows;
		this->m_numRows	= this->m_numCols;
		this->m_numCols	= r;

		for( int i=0; i<this->m_numRows; ++i )
		{
			for( int j=0; j<this->m_numCols; ++j )
			{
				(*this)(i, j) = *(temp + this->m_numRows*j + i);
			}// end of j loop
		}// end of i loop

		SafeDelete( temp );
	}


	// https://softwareengineering.stackexchange.com/questions/271713/transpose-a-matrix-without-a-buffering-one
	/// Performs in-place transposition of this.
	void TransposeInplace()
	{
		int m = this->m_numRows;
		int n = this->m_numCols;

		for( int start=0; start<this->numMaxElms(); ++start )
		{
			int next = start;
			int i = 0;
			do
			{
				++i;
				next = (next % m) * n + next / m;
			} while (next > start);

			if (next >= start && i != 1)
			{	
				const T tmp = (*this)(start);//const double tmp = m[start];
				next = start;
				do
				{
					i = (next % m) * n + next / m;
					(*this)(next) = i==start ? tmp : (*this)(i);//m[next] = (i == start) ? tmp : m[i];
					next = i;
				} while (next > start);
			}
		}

		this->m_numRows = n;
		this->m_numCols = m;

	}




	//=========================== virtual methods =============================//

	T *ptr( int row, int col ) const// row: vertical, col: horizontal
	{
		return (T *)( m_pData + this->m_numCols*row + col );
	}


	T *ptr( int i ) const
	{
		return (T *)( m_pData + i );
	}


	unsigned char *pData() const
	{
		return (unsigned char *)m_pData;
	}



protected:

	T	*m_pData;	//OreOreLib::Memory<T>	m_Data;



	void Allocate( int numrows, int numcols )
	{
		assert( numrows>0 && numcols>0 );

		Release();

		this->m_numRows	= numrows;
		this->m_numCols	= numcols;
		m_pData			= new T[ numrows * numcols ]();
	}


};




#endif	// DYNAMIC_MATRIX_H //