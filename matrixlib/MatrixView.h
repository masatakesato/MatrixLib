#ifndef MATRIX_VIEW_H
#define	MATRIX_VIEW_H

#include	<oreore/mathlib/MathLib.h>

#include	"IMatrix.h"



//##########################################################################//
//							Matrix View Class								//
//##########################################################################//

template< typename T >	struct VIEW{ using Type = typename T; };


template< typename >
struct is_view : std::false_type{};

template< typename T >
struct is_view< VIEW<T> > : std::true_type{};

template< typename T >
constexpr bool is_view_v = is_view<T>::value;




template< typename T >
using MatrixView = Matrix< VIEW<T>, DynamicDim, DynamicDim >;



template< typename T >
using RowVectorView = Matrix< VIEW<T>, 1, DynamicDim >;


template< typename T >
using ColVectorView = Matrix< VIEW<T>, DynamicDim, 1 >;




template< typename V, int R, int C >
class Matrix< V, R, C, std::enable_if_t< is_view_v<V> && (R==DynamicDim || C==DynamicDim) > > : public IMatrix< typename V::Type >
{
	using T = typename V::Type;
public:

	// Default constructor
	Matrix()
		: IMatrix<T>( MEM::PARTIAL )
		, m_refData( nullptr )
		, m_ColOffset( 0 )
	{
	}


	// Constructor using IMatrix
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	Matrix( const IMatrix<T>& mat, int i, int j, int numrows, int numcols )
		: IMatrix<T>( numrows, numcols, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.numCols() )
	{
	}


	// Constructor using MatrixView
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	Matrix( const Matrix& mat, int i, int j, int numrows, int numcols )
		: IMatrix<T>( numrows, numcols, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.m_ColOffset )
	{
	}


	// Constructor using array
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	Matrix( const T pdata[], int numrows, int numcols )
		: IMatrix<T>( numrows, numcols, MEM::PARTIAL )
		, m_refData( pdata )
		, m_ColOffset( 0 )
	{
		//assert( numrows * numcols > pdata.Length() );
	}


	// Copy constructor
	Matrix( const Matrix& obj )
		: IMatrix<T>( obj.m_numRows, obj.m_numCols, MEM::PARTIAL )
		, m_refData( obj.m_refData )
		, m_ColOffset( obj.m_ColOffset )
	{
	}


	// Move constructor
	Matrix( Matrix&& obj )
		: IMatrix<T>( obj.m_numRows, obj.m_numCols, MEM::PARTIAL )
		, m_refData( obj.m_refData )
		, m_ColOffset( obj.m_ColOffset )
	{
		obj.m_numRows	= 0;
		obj.m_numCols	= 0;	
		obj.m_refData	= nullptr;// clear reference from obj
		obj.m_ColOffset	= 0;
	}


	// Destructor
	~Matrix()
	{
		m_refData = nullptr;
	}

	
	// Copy assignment operator
	Matrix& operator=( const Matrix& obj )
	{
		//tcout << _T("//========= MatrixView Copy assignment operator... ==========//\n" );
		if( this != &obj )
		{
			this->m_numRows = obj.m_numRows;
			this->m_numCols = obj.m_numCols;
			m_refData		= obj.m_refData;
			m_ColOffset		= obj.m_ColOffset;
		}

		return *this;
	}


	// Move assignment operator for MatrixView
	Matrix& operator=( Matrix&& obj )
	{
		if( this != &obj )
		{
			this->m_numRows	= obj.m_numRows;
			this->m_numCols	= obj.m_numCols;
			m_refData		= obj.m_refData;
			m_ColOffset		= obj.m_ColOffset;

			// clear reference from src
			obj.m_numRows	= 0;
			obj.m_numCols	= 0;
			obj.m_refData	= nullptr;
			obj.m_ColOffset	= 0;
		}

		return *this;
	}


	// MatrixView-specific implementation (RHS must be square matrix)
	template< typename U >
	std::enable_if_t<std::is_arithmetic_v<U>, Matrix& >
	operator*=( const IMatrix<U>& right )
	{
		assert( (this->m_numCols==right.numRows()) && (this->m_numCols==right.numCols()) );

		DynamicMatrix<U> out( this->m_numRows, this->m_numCols );

		for( int i=0; i<out.numRows(); ++i )
		{
			for( int j=0; j<out.numCols(); ++j )
			{
				for( int k=0; k<this->m_numCols; ++k )
					out(i, j) += (*this)(i, k) * right(k, j);
			}
		}

		this->Copy( out );

		return *this;
		//return *this = (*this) * right;// MatrixViewだけじゃムリ
	}


	// 禁止. 新規作成する戻り値は領域持ってないとダメ(Matrix or StaticMatrixのみ可 )
	//MatrixView<T> operator+( const MatrixView& right )
	//{
	//	assert( right.IsSameSize(*this) );

	//	MatrixView<T> m(*this);// copy constructor of Matrix class
	//	//return static_cast<MatrixView<T>&>( m);
	//	return  static_cast<MatrixView<T>&>( m += static_cast<MatrixView<T>>(right) );
	//}


	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	void Init( const Matrix& mat, int i, int j, int numrows, int numcols )
	{
		this->m_numRows = numrows;
		this->m_numCols = numcols;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.m_ColOffset;
	}


	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	void Init( const IMatrix<T>& mat, int i, int j, int numrows, int numcols )
	{
		this->m_numRows = numrows;
		this->m_numCols = numcols;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.numCols();
	}


	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	void Init( const T& pdata, int numrows, int numcols )
	{
		this->m_numRows = numrows;
		this->m_numCols = numcols;
		m_refData		= pdata;
		m_ColOffset		= 0;
	}


	void Release()
	{
		this->m_numRows	= 0;
		this->m_numCols	= 0;
		m_refData		= nullptr;
		m_ColOffset		= 0;
	}


	void Clear()
	{
		for( int i=0; i<this->m_numRows; ++i )
			memset( (T*)m_refData + m_ColOffset*i, 0, this->m_numCols * sizeof(T) );
	}


	bool Copy( const IMatrix<T>& src )
	{
		if( !src.IsSameSize(*this) )
			return false;

		for( int i=0; i<this->m_numRows; ++i )
			memcpy_s( ptr(i,0), this->m_numCols * sizeof(T), src.ptr(i,0), src.numCols() * sizeof(T) );

		return true;
	}


	bool CopyTo( IMatrix<T>& dst ) const
	{
		if( !dst.IsSameSize(*this) )
			return false;

		for( int i=0; i<dst.numRows(); ++i )
			memcpy_s( dst.ptr(i,0), dst.numCols() * sizeof(T), ptr(i,0), this->m_numCols * sizeof(T) );

		return true;
	}


	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==DynamicDim>* = nullptr >
	void Transpose()
	{
		assert( this->IsSquare() );

		T tmp;

		for( int i=0; i<this->m_numRows; ++i )
		{
			for( int j=i+1; j<this->m_numCols; ++j )
			{
				tmp = *ptr(i, j);
				*ptr(i, j) = *ptr(j, i);
				*ptr(j, i) = tmp;
			}// end of j loop
		}// end of i loop

	}



	//=========================== virtual methods =============================//

	T *ptr( int row, int col ) const// row: vertical, col: horizontal
	{
		return (T *)( m_refData + m_ColOffset * row + col );
	}


	T *ptr( int i ) const
	{
		static int div, mod;
		DivMod( div, mod, i, this->m_numCols );
		return (T *)(m_refData + m_ColOffset * div + mod);
	}


	unsigned char *pData() const
	{
		return (unsigned char *)m_refData;
	}



	//===================== RowVectorView specific methods ====================//

	// Constructor using IMatrix
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_==DynamicDim>* = nullptr >
	Matrix( const IMatrix<T>& mat, int i, int j, int numcols )
		: IMatrix<T>( R_, numcols, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.numCols() )
	{
	}

	// Constructor using MatrixView
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_==DynamicDim>* = nullptr >
	Matrix( const Matrix& mat, int i, int j, int numcols )
		: IMatrix<T>( R_, numcols, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.m_ColOffset )
	{
	}

	// Constructor using array
	template< int R_=R, int C_=C, std::enable_if_t<R_==1 && C_==DynamicDim>* = nullptr >
	Matrix( const T pdata[], int numcols )
		: IMatrix<T>( R_, numcols, MEM::PARTIAL )
		, m_refData( pdata )
		, m_ColOffset( 0 )
	{
	}


	// Init using IMatrix
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==1 && C_==DynamicDim, void >
	Init( const IMatrix<T>& mat, int i, int j, int numcols )
	{
		this->m_numRows = R_;
		this->m_numCols = numcols;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.numCols();
	}

	// Init using MatrixView
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==1 && C_==DynamicDim, void >
	Init( const Matrix& mat, int i, int j, int numcols )
	{
		this->m_numRows = R_;
		this->m_numCols = numcols;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.m_ColOffset;
	}

	// Init using array
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==1 && C_==DynamicDim, void >
	Init( const T& pdata, int numcols )
	{
		this->m_numRows = R_;
		this->m_numCols = numcols;
		m_refData		= pdata;
		m_ColOffset		= 0;
	}



	//===================== ColVectorView specific methods ====================//
	
	// Constructor using IMatrix
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==1>* = nullptr >
	Matrix( const IMatrix<T>& mat, int i, int j, int numrows )
		: IMatrix<T>( numrows, C_, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.numCols() )
	{
	}

	// Constructor using MatrixView
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==1>* = nullptr >
	Matrix( const Matrix& mat, int i, int j, int numrows )
		: IMatrix<T>( numrows, C_, MEM::PARTIAL )
		, m_refData( mat.ptr(i, j) )
		, m_ColOffset( mat.m_ColOffset )
	{
	}

	// Constructor using array
	template< int R_=R, int C_=C, std::enable_if_t<R_==DynamicDim && C_==1>* = nullptr >
	Matrix( const T pdata[], int numrows )
		: IMatrix<T>( numrows, C_, MEM::PARTIAL )
		, m_refData( pdata )
		, m_ColOffset( 0 )
	{
	}


	// Init using IMatrix
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==DynamicDim && C_==1, void >
	Init( const IMatrix<T>& mat, int i, int j, int numrows )
	{
		this->m_numRows = numrows;
		this->m_numCols = C_;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.numCols();
	}

	// Init using MatrixView
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==DynamicDim && C_==1, void >
	Init( const Matrix& mat, int i, int j, int numrows )
	{
		this->m_numRows = numrows;
		this->m_numCols = C_;
		m_refData		= mat.ptr(i, j);
		m_ColOffset		= mat.m_ColOffset;
	}

	// Init using array
	template< int R_=R, int C_=C >
	std::enable_if_t< R_==DynamicDim && C_==1, void >
	Init( const T& pdata, int numrows )
	{
		this->m_numRows = numrows;
		this->m_numCols = C_;
		m_refData		= pdata;
		m_ColOffset		= 0;
	}




protected:

	const T*	m_refData;
	int			m_ColOffset;

};



#endif	// MATRIX_VIEW_H //