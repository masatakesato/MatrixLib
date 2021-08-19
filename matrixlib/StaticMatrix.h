#ifndef STATIC_MATRIX_H
#define	STATIC_MATRIX_H

#include	<oreore/meta/TypeTraits.h>

#include	"IMatrix.h"




//##########################################################################//
//						N x M Static Matrix Class							//
//##########################################################################//


template< typename T, size_t R, size_t C >
using StaticMatrix = Matrix< T, R, C, std::enable_if_t< R!=DynamicDim && C!=DynamicDim > >;



template< typename T, size_t C >
using StaticRowVector = Matrix< T, 1, C, std::enable_if_t< C!=DynamicDim > >;



template< typename T, size_t R >
using StaticColVector = Matrix< T, R, 1, std::enable_if_t< R!=DynamicDim > >;





template< typename T, size_t R, size_t C >
class Matrix< T, R, C, std::enable_if_t< std::is_arithmetic_v<T> && R!=DynamicDim && C!=DynamicDim > > : public IMatrix<T>
{
public:

	// Default constructor
	Matrix()
		: IMatrix<T>(R, C, MEM::ENTIRE )
	{
		this->Clear();
	}


	// Constructor
	Matrix( const T* const data )
		: IMatrix<T>( R, C, MEM::ENTIRE )
	{
		//assert( numrow>0 && numcol>0 );
		memcpy_s( m_Data, this->Size(), data, this->Size() );
	}

	 
	// Constructor
	template< int R_=R, int C_=C, std::enable_if_t< R_!=DynamicDim && C_!=DynamicDim >* = nullptr >
	Matrix( std::initializer_list<std::initializer_list<T>> ilist )
		: IMatrix<T>( R_, C_, MEM::ENTIRE )
	{
		// Extract numRows and numCols from ilist
		assert( R_ == int( ilist.size() ) );
		for( auto &row_data : ilist )
			assert( C_ == int( row_data.size() ) );

		// copy value
		auto p = &m_Data[0];
		for( auto& row_data : ilist )
			for( const auto& val : row_data )
				*p++ = val;
	}


	// Constructor for row vector
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_!=DynamicDim>* = nullptr >
	Matrix( std::initializer_list<T> ilist )
		: IMatrix<T>( R_, C_, MEM::ENTIRE )
	{
		assert( int(ilist.size()) == C_ );

		auto p = &m_Data[0];
		for( const auto& val : ilist )	*p++ = val;
	}


	// Constructor for col vector
	template< int R_=R, int C_=C, std::enable_if_t<R_!=DynamicDim && C_==1>* = nullptr >
	Matrix( std::initializer_list<T> ilist )
		: IMatrix<T>( R_, C_, MEM::ENTIRE )
	{
		assert( int(ilist.size()) == R_ );

		auto p = &m_Data[0];
		for( const auto& val : ilist )	*p++ = val;
	}


	// Constructor
	Matrix( const IMatrix<T> &obj )
		: IMatrix<T>( R, C, MEM::ENTIRE )
	{
		if( obj.IsSameSize(*this) )
			Copy( obj );
		else
			this->Clear();
	}


	// Copy constructor
	Matrix( const Matrix<T,R,C>& obj )
		: IMatrix<T>( R, C, MEM::ENTIRE )
	{
		memcpy_s( m_Data, this->Size(), obj.pData(), obj.Size() );
	}



	~Matrix()
	{
		//this->m_pData = nullptr;
	}


	// Copy assignment operator
	Matrix& operator=( const IMatrix<T>& obj )
	{
		assert( obj.IsSameSize(*this) );

		if( this != &obj && obj.IsSameSize(*this) )
			Copy( obj );//memcpy( m_Data, obj.pData(), (size_t)(size_t)(R * C * sizeof(T)) );

		return *this;
	}



	// Static Matrix-specific implementation (RHS must be square matrix)
	template< typename U >
	std::enable_if_t<std::is_arithmetic_v<U>, Matrix& >
	operator*=( const IMatrix<U>& right )
	{
		assert( (this->m_numCols==right.numRows()) && (this->m_numCols==right.numCols()) );

		static Matrix<U,R,C> out;

		for( int i=0; i<out.m_numRows; ++i )
		{
			for( int j=0; j<out.m_numCols; ++j )
			{
				out(i, j) = 0;
				for( int k=0; k<this->m_numCols; ++k )
					out(i, j) += (*this)(i, k) * right(k, j);
			}
		}

		this->Copy( out );

		return *this;
		//return *this = (*this) * right;// MatrixViewだけじゃムリ
	}


	void Clear()
	{
		memset( m_Data, 0, this->Size() );
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


	void Transpose()
	{
		T *temp	= new T[ this->numElms() ];
		memcpy_s( temp, this->Size(), m_Data, this->Size() );

		int r			= this->m_numRows;
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
		return (T *)( m_Data + this->m_numCols*row + col );
	}


	T *ptr( int i ) const
	{
		return (T *)( m_Data + i );
	}


	unsigned char *pData() const
	{
		return (unsigned char *)m_Data;
	}



protected:

	T	m_Data[ R * C ];
};




#endif	// STATIC_MATRIX_H //