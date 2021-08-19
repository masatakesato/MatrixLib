#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#include <matrixlib/MatrixLib.h>
#include <matrixlib/GaussElimination.h>
//#include <matrixlib/ICA.h>



int main(int argc, char *argv[])
{
/*
	srand(time(NULL));
	// y = Wx
	Matrix<double> y( 5, 2 );// 独立成分ベクトル(Wとxから推定)
	Matrix<double> W( y.numRows(), y.numRows() );// なんだか分からない重み行列
	Matrix<double> x( 5, 2 );// ベクトル(計測データ)
	// x1 =Ay
	Matrix<double> A( y.numRows(), y.numRows() );
	Matrix<double> x1( x.numRows(), x.numCols() );// ベクトル(復元計測データ)

	for( int i=0; i<x.numRows(); i++ )
	{
		for( int j=0; j<x.numCols(); j++ )
		{
			x(i,j) = double(rand()%255) / 255.0;
		}
	}

	//for( int i=0; i<x.numElms(); i++ )
	//{
	//	x()
	//	*(x->Elements + i) = rand()%255;
	//	*(x->Elements + i) /= 255;
	//	//*(x->Elements + i) *= 2;
	//	//*(x->Elements + i) -= 1;
	//}

	//ICA_Jutten91(y,W,x);
	//
	//Gauss_Inverse(A,W);//Wの逆行列をAに格納
	//Matrix_mult(x1,A,y);//独立成分yと重み行列Aを使って計測データを復元

	//printf("W:\n");
	//Matrix_display(W);
	//
	//printf("x:\n");
	//Matrix_display(x);
	//
	//printf("x1:\n");
	//Matrix_display(x1);
*/
return 0;
}