#include <matrixlib/MatrixLib.h>
#include <matrixlib/Covariance.h>




int main( int argc, char *argv[] )
{
	OreOreLib::init_genrand( 0 );

	DynamicMatrix<double>	a( 5, 3 );							// 3 x 2 matrix 
	DynamicMatrix<double>	a_cov( a.numCols(), a.numCols() );	// covariance mateix of a

	//Random( a, 0.0, 1.0 );// Initialize a with [0, 1] random variables.

	a(0, 0) = 90;	a(0, 1) = 60;	a(0, 2) = 90;
	a(1, 0) = 90;	a(1, 1) = 90;	a(1, 2) = 30;
	a(2, 0) = 60;	a(2, 1) = 60;	a(2, 2) = 60;
	a(3, 0) = 60;	a(3, 1) = 60;	a(3, 2) = 90;
	a(4, 0) = 30;	a(4, 1) = 30;	a(4, 2) = 30;

//	a.Transpose();
	a.Display();

	
	// Calc covariance matrix
	Covariance( a_cov, a );// aの分散共分散行列をa_covに格納

	a_cov.Display();

//	printf("Matrix a:\n");
//	Matrix_display(a);
//	printf("\nVariance-Covariance Matrix:\n");
//	Matrix_display(a_cov);

	return 0;
}