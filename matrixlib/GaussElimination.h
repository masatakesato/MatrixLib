#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H


#include	"MatrixLib.h"



template< typename T >
void Gauss_Elimination( IMatrix<T>& inout )
{

	//============== Forward elimination: make lower triangulat part to zero ============//
	for( int j=0; j<inout.numCols(); j++ )// 列単位でゼロを充填していく
	{
		// Pivot selection
		int row_selected = j;
		for( int k=j; k<inout.numRows(); k++ )
		{
			if( fabs( *inout.ptr(k, j) ) > fabs( *inout.ptr(row_selected, j) ) )
				row_selected = k;
		}
		if( row_selected != j )// j行目とrow_selected行目を入れ替える
		{
			for( int l=0; l<inout.numCols(); l++ )
			{
				T tmp = *inout.ptr( j, l );
				*inout.ptr(j, l) = *inout.ptr(row_selected, l);
				*inout.ptr(row_selected, l) = tmp;
			}
		}

		for( int i=j+1; i<inout.numRows(); i++ )// 各行について，，，
		{
			// i行目からj行ベクトルを減算する
			T scale = *inout.ptr(i, j) / *inout.ptr(j, j);
			for( int l=0; l<inout.numCols(); l++ )
				*inout.ptr(i, l) -= scale * *inout.ptr(j, l);

		}// end of i loop


	}// end of j loop


	 //=============== 後退代入：行列の最下行から順次変数の値を決めていく ===============//
	for( int i=inout.numRows()-1; i>=0; i-- )
	{
		// 自分より上の行ベクトルに対して減算処理を行う
		for( int k=i-1; k>=0; k-- )
		{
			T scale = *inout.ptr(k, i) / *inout.ptr(i, i);
			for( int j=0; j<inout.numCols(); j++ )
				*inout.ptr(k, j) -= scale * *inout.ptr(i, j);
		}

		// i行目の行ベクトルを正規化する
		T scale = 1.0 / *inout.ptr(i, i);
		for( int j=0; j<inout.numCols(); j++ )
			*inout.ptr(i, j) *= scale;

	}// end of i loop

}



template< typename T >
void Inverse_Gauss( IMatrix<T>& out, const IMatrix<T>& in )
{
	assert( in.numRows()==out.numRows() && in.numCols()==out.numCols() );

	DynamicMatrix<T> buf( in.numRows(), in.numCols()*2 ); // buf = [ in | I ]

	for( int i=0; i<in.numRows(); ++i )
	{
		// fill original matrix elements on left part
		for( int j=0; j<in.numCols(); ++j )
			buf(i, j) = in(i, j);

		// set idendity matrix on right part
		buf(i, i+in.numCols()) = (T)1;
	}

	// execute gaussian elimination
	Gauss_Elimination(buf);

	// copy result to out matrix
	for( int i=0; i<out.numRows(); ++i )
		for( int j=0; j<out.numCols(); ++j )
			out(i, j) = buf(i, j+out.numCols());

}




#endif /* GAUSS_ELIMINATION_H */