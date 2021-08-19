#ifndef MATRIX_OPERATIONS_H
#define	MATRIX_OPERATIONS_H

#include	<oreore/mathlib/MathLib.h>
#include	<oreore/mathlib/MersenneTwister.h>
#include	<oreore/container/Array.h>



//######################################################################################//
//																						//
//									Initialization										//
//																						//
//######################################################################################//


// Clear
template< typename T >
inline void Zero( IMatrix<T>& inout )
{
	inout.Clear();
}



// Fill all elements with one
template< typename T >
inline void One( IMatrix<T>& inout )
{
	for( int i=0; i<inout.numElms(); ++i )
		inout(i) = (T)1;
}



// Idenity
template< typename T >
inline void Identity( IMatrix<T>& inout )
{
	inout.Clear();
	for( int i=0; i<Min(inout.numRows(), inout.numCols()); ++i )
		inout(i, i) = (T)1;
}



// Randomize
template< typename T >
inline void Random( IMatrix<T>& inout, T min=0, T max=1 )
{
	for(int row=0; row<inout.numRows(); ++row)
	{
		for(int col=0; col<inout.numCols(); ++col)
		{
			inout( row, col ) = OreOreLib::Rand<T>( min, max );
		}// end of col loop
	}// end of row loop
}



// Transpose
template< typename T >
inline void Transpose( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( in.numRows()!=out.numCols() || in.numCols()!=out.numRows() )
		return;

	for( int i=0; i<out.numRows(); ++i )
	{
		for(int j=0; j<out.numCols(); ++j )
		{
			out(i, j) = in(j, i);
		}
	}
}



template< typename T >
inline void InitRow( IMatrix<T>& inout, const T data[], int row=0 )
{
	for( int i=0; i<inout.numCols(); ++i )
	{
		inout(row, i) = data[i];
	}



}



/*
// out = in^-1
template< typename T >
inline void Inverse( DynamicMatrix<T>& out, const DynamicMatrix<T>& in ) 
{	
	if( !in.IsSquare() )
		return;

	T buf; //一時的なデータを蓄える
	int i,j,k; //カウンタ
	int	N = in.numRows();  //配列の次数	

						   // 初期化
	if( out.IsEmpty() ) out.Init( in.numRows(), in.numCols() );
	Identity( out );

	//掃き出し法
	for(i=0; i<N; i++)
	{
		buf = in(i, i)>1.0e-8 ? 1.0 / in(i, i) : 0.0;
		buf = buf > 1.0e-8 ? buf : 0.0;

		for(j=0; j<N; j++)
		{
			in(i, j) *= buf;
			out(i, j) *= buf;
		}

		for(j=0; j<N; j++)
		{
			if(i!=j)
			{
				buf = in(j, i);
				for(k=0; k<N; k++)
				{
					in(j, k) -= in(i, k) * buf;
					out(j, k) -= out(i, k) * buf;
				}
			}
		}// end of j loop
	}// end of i loop
}



// Moore Penrose pseudo inverse
template< typename T >
void PresudoInverse( DynamicMatrix<T> &A_pinv, const DynamicMatrix<T> &A )
{

	DynamicMatrix<T>	At;
	DynamicMatrix<T>	AtA;

	Transpose( At, A );

	Multiply( AtA, At, A );

	//AtA.Display();

	DynamicMatrix<T>	AtA_inv;
	Inverse( AtA_inv, AtA );

	AtA_inv.Display();

	Multiply( A_pinv, AtA_inv, At );
}
*/



//######################################################################################//
//																						//
//								Arithmetic Fuinctions									//
//																						//
//######################################################################################//


// out = in1 + in2
template< typename T >
inline void Add( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( !out.IsSameSize(in1) || !out.IsSameSize(in2) ) return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) = in1(i) + in2(i);
}



// out = in1 + scalar
template< typename T >
inline void Add( IMatrix<T>& out, const IMatrix<T>& in, const T& scalar )
{
	for( int i=0; i<out.numElms(); ++i )
		out(i) = in(i) + scalar;
}



// out += in
template< typename T >
inline void Add( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( !out.IsSameSize(in) )	return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) += in(i);
}



// out += scalar
template< typename T >
inline void Add( IMatrix<T>& inout, const T& scalar )
{
	for( int i=0; i<inout.numElms(); ++i )
		inout(i) += scalar;
}



// out = in1 * scalar1 + in2 * scalar2
template< typename T >
inline void AddScaled( IMatrix<T>& out, const IMatrix<T>& in1, const T& scalar1, const IMatrix<T>& in2, const T& scalar2 )
{
	if( !out.IsSameSize(in1) || !out.IsSameSize(in2) ) return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) = in1(i) * scalar1 + in2(i) * scalar2;
}



// out += in * scalar
template< typename T >
inline void AddScaled( IMatrix<T>& out, const IMatrix<T>& in, const T& scalar )
{
	if( !out.IsSameSize(in) )	return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) += in(i) * scalar;
}



// out = in1 - in2
template< typename T >
inline void Subtract( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( !out.IsSameSize(in1) || !out.IsSameSize(in2) ) return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) = in1(i) - in2(i);
}



// out -= in
template< typename T >
inline void Subtract( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( !out.IsSameSize(in) )	return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) -= in(i);
}



// out = in1 * in2
template< typename T >
inline void Multiply( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( out.numRows()!=in1.numRows() || out.numCols()!=in2.numCols() || in1.numCols()!=in2.numRows() ) return;

	for( int row=0; row<out.numRows(); ++row )
	{
		for( int col=0; col<out.numCols(); ++col )
		{
			T& dest	= out(row, col);
			dest = 0;
			for( int k=0; k<in1.numCols(); ++k )
				dest += in1(row, k) * in2(k, col);
		}
	}

}



// DynamicMatrix<T> *= IMatrix<T>
template< typename T >
inline void Multiply( DynamicMatrix<T>& out, const IMatrix<T>& in )
{
	out *= in;
}


// Static Matrix *= any type of matrix
template< typename T, size_t R, size_t C >
inline void Multiply( /*Static*/Matrix<T, R, C>& out, const IMatrix<T>& in )
{
	out *= in;
}


// MatrixView *= any type of matrix
template< typename T >
inline void Multiply( MatrixView<T>/*Matrix<VIEW<T>>*//*View<T>*/& out, const IMatrix<T>& in )
{
	out *= in;
}



// out = scalar * in1
template< typename T >
inline void Scale( IMatrix<T>& out, const T& scalar, const IMatrix<T>& in1 )
{
	if( !in1.IsSameSize(out) )	return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) = scalar * in1(i);
}



// inout *= scalar
template< typename T >
inline void Scale( IMatrix<T>& inout, const T& scalar )
{
	for( int i=0; i<inout.numElms(); ++i )
		inout(i) *= scalar;
}



// out = in1 * in2^T
template< typename T >
inline void MultiplyTranspose( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( in1.numCols()!=in2.numCols() || out.numRows()!=in1.numRows() || out.numCols()!=in2.numRows() )
		return;

	for( int row=0; row<out.numRows(); ++row )
	{
		for( int col=0; col<out.numCols(); ++col )
		{
			T& dest	= out(row, col);
			dest	= 0;

			for( int k=0; k<in1.numCols(); ++k )
				dest += in1(row, k) * in2(col, k);

		}// end of col loop
	}// end of row loop

}



// out = in * in^T
template< typename T >
inline void MultiplyTranspose( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( out.numRows()!=in.numRows() || out.numCols()!=in.numRows() )
		return;

	for( int row=0; row<out.numRows(); ++row )
	{
		for( int col=0; col<out.numCols(); ++col )
		{
			T& dest	= out(row, col);
			dest	= 0;

			for( int k=0; k<in.numCols(); ++k )
				dest += in(row, k) * in(col, k);

		}// end of col loop
	}// end of row loop

}



// out = in1^T *  in2
template< typename T >
inline void TransposeMultiply( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( in1.numRows()!=in2.numRows() || out.numRows()!=in1.numCols() || out.numCols()!=in2.numCols() )
		return;

	for( int row=0; row<out.numRows(); ++row )
	{
		for( int col=0; col<out.numCols(); ++col )
		{
			T& dest	= out(row, col);
			dest	= 0;

			for( int k=0; k<in1.numRows(); ++k )
				dest += in1(k, row) * in2(k, col);

		}// end of col loop
	}// end of row loop

}



// out = in^T *  in
template< typename T >
inline void TransposeMultiply( IMatrix<T>& out, const IMatrix<T>& in )
{
	if( out.numRows()!=in.numCols() || out.numCols()!=in.numCols() )
		return;

	for( int row=0; row<out.numRows(); ++row )
	{
		for( int col=0; col<out.numCols(); ++col )
		{
			T& dest	= out(row, col);
			dest	= 0;

			for( int k=0; k<in.numRows(); ++k )
				dest += in(k, row) * in(k, col);

		}// end of col loop
	}// end of row loop

}



// out = Hadamard product of in1 and in2
template< typename T >
inline void Hadamard( IMatrix<T>& out, const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	if( !out.IsSameSize(in1) || !in1.IsSameSize(in2) )
		return;

	for( int i=0; i<out.numElms(); ++i )
		out(i) = in1(i) * in2(i);
}



// Trace(sum of diagonal elements) of matrix
template< typename T >
inline T Tr( const IMatrix<T>& in )
{
	if( !in.IsSquare() )
		return nan("");

	T val = 0;
	for( int i=0; i<in.numRows(); ++i )
		val += in(i, i);

	return val;
}




//######################################################################################//
//																						//
//									Partial Calculation									//
//																						//
//######################################################################################//


// Dot product. ( in1/in2 must be (1xn) or (n*1) )
template< typename T >
inline T Dot( const IMatrix<T>& in1, const IMatrix<T>& in2 )
{
	T dot_product = 0;
	for( int i=0; i<in1.numElms(); ++i )
		dot_product += in1(i) * in2(i);

	return dot_product;
}



// Dot product of row by column ( in1[row1, 1:] and in2[1:, col2] ).
template< typename T >
inline T DotRC( const IMatrix<T>& in1, const IMatrix<T>& in2, int row1=0, int col2=0 )
{
	assert( in1.numCols()==in2.numRows() && row1>=0 && row1<in1.numRows() && col2>=0 && col2<in2.numCols() );

	T dot_product = 0;
	for( int k=0; k<in1.numCols(); ++k )
		dot_product += in1(row1, k) * in2(k, col2);

	return dot_product;
}



// Dot product of row by row ( in1[1:, col1] and in2[1:, col2] ).
template< typename T >
inline T DotRR( const IMatrix<T>& in1, const IMatrix<T>& in2, int row1=0, int row2=0 )
{
	assert( in1.numCols()==in2.numCols() && row1>=0 && row1<in1.numRows() && row2>=0 && row2<in2.numRows() );

	int offset1 =row1 * in1.numCols();
	int offset2 =row2 * in2.numCols();

	T dot_product = 0;
	for( int i=0; i<in1.numCols(); ++i )
		dot_product += in1(offset1 + i) * in2(offset2 + i);

	return dot_product;
}



// Dot product of column by column ( in1[1:, col1] and in2[1:, col2] ).
template< typename T >
inline T DotCC( const IMatrix<T>& in1, const IMatrix<T>& in2, int col1=0, int col2=0 )
{
	assert( in1.numRows()==in2.numRows() && col1>=0 && col1<in1.numCols() && col2>=0 && col2<in2.numCols() );

	int offset = in1.numCols();

	T dot_product = 0;
	for( int i=0; i<in1.numRows(); ++i )
		dot_product += in1(col1 + i*offset) * in2(col2 + i*offset);

	return dot_product;
}



// Sum of single matrix row
template< typename T >
inline T SumRow( const IMatrix<T>& in, int row=0 )
{
	assert( row>=0 && row<in.numRows() );

	int offset = in.numCols() * row;

	T sum = 0;
	for( int i=0; i<in.numCols(); ++i )
		sum += in(offset + i);

	return sum;
}



// Sum of single matrix column
template< typename T >
inline T SumCol( const IMatrix<T>& in, int col=0 )
{
	assert( col>=0 && col<in.numCols() );

	int offset = in.numCols();

	T sum = 0;
	for( int i=0; i<in.numRows(); ++i )
		sum += in(col + i*offset);

	return sum;
}




//######################################################################################//
//																						//
//										Permutation										//
//																						//
//######################################################################################//

// Permutation
template< typename T >
void Permutation( IMatrix<T>& inout, int num, const int perm[] )
{
	assert( inout.IsSquare() && inout.numRows()==num );
	Zero( inout );
	for( int i=0; i<inout.numRows(); ++i )	inout(i, perm[i]) = (T)1;
}


template< typename T >
void Permutation( IMatrix<T>& inout, const OreOreLib::Array<int>& perm )
{
	assert( inout.IsSquare() && inout.numRows()==perm.Length() );
	Zero( inout );
	for( int i=0; i<inout.numRows(); ++i )	inout(i, perm[i]) = (T)1;
}



// Transposed Permutation
template< typename T >
void TransposedPermutation( IMatrix<T>& inout, int num, const int perm[] )
{
	assert( inout.IsSquare() && inout.numCols()==num );
	Zero( inout );
	for( int i=0; i<inout.numCols(); ++i )	inout(perm[i], i) = (T)1;
}


template< typename T >
void TransposedPermutation( IMatrix<T>& inout, const OreOreLib::Array<int>& perm )
{
	assert( inout.IsSquare() && inout.numCols()==perm.Length() );
	Zero( inout );
	for( int i=0; i<inout.numCols(); ++i )	inout(perm[i], i) = (T)1;
}





//######################################################################################//
//																						//
//								Arithmetic Operators									//
//																						//
//######################################################################################//


template< typename T >
const DynamicMatrix<T> operator*( const IMatrix<T>& left, const IMatrix<T>& right )
{
	assert( left.numCols()==right.numRows() );

	DynamicMatrix<T> out( left.numRows(), right.numCols() );

	for( int i=0; i<out.numRows(); ++i )
	{
		for( int j=0; j<out.numCols(); ++j )
		{
			T& dest	= out(i, j);
			dest = 0;
			for(int k=0; k<left.numCols(); ++k)
				dest += left(i, k) * right(k, j);
		}
	}

	return out;
}



template< typename T >
const DynamicMatrix<T> operator*( const T& scalar, const IMatrix<T>& right ) 
{ 
	DynamicMatrix<T> out( right.numRows(), right.numCols() );

	for( int i=0; i<out.numElms(); ++i )
		out(i) = scalar * right(i);

	return out;
} 



template< typename T >
const DynamicMatrix<T> operator+( const IMatrix<T>& left, const IMatrix<T>& right )
{
	DynamicMatrix<T> m(left.numRows(), left.numCols());

	for( int i=0; i<m.numElms(); ++i )
		m(i) = left(i) + right(i);

	return m;
}



template< typename T >
const DynamicMatrix<T> operator-( const IMatrix<T>& left, const IMatrix<T>& right )
{
	DynamicMatrix<T> m(left.numRows(), left.numCols());

	for( int i=0; i<m.numElms(); ++i )
		m(i) = left(i) - right(i);

	return m;
}




#endif // !MATRIX_OPERATIONS_H
