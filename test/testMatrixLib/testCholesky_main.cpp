#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/Cholesky.h>






int main( int argc, char *argv[] )
{
/*
// https://en.wikipedia.org/wiki/Cholesky_decomposition
	Matrix<double> A( 3, 3 );	// Matrix to be decomposed.
	Matrix<double> L( 3, 3 );	// Lower triangular matrix
	Matrix<double> Lt( 3, 3 );	// Upper triangular matrix

	// Init Matrix A
	A(0, 0) = 4;		A(0, 1) = 12;	A(0, 2) = -16;
	A(1, 0) = 12;	A(1, 1) = 37;	A(1, 2) = -43;
	A(2, 0) = -16;	A(2, 1) = -43;	A(2, 2) = 98;
	//Matrix_insert_at(0,0,4,A); Matrix_insert_at(0,1,2,A); Matrix_insert_at(0,2,2,A);
	//Matrix_insert_at(1,0,2,A); Matrix_insert_at(1,1,5,A); Matrix_insert_at(1,2,3,A);
	//Matrix_insert_at(2,0,2,A); Matrix_insert_at(2,1,3,A); Matrix_insert_at(2,2,11,A);

	A.Display();

	// Decompose A to L, Lt
	Cholesky( L, A );
	Transpose( Lt, L );

	tcout << "L:\n";
	L.Display();

	tcout << "Lt\n";
	Lt.Display();


//	Matrix<double> U( 3, 3 );	// U = L*Lt

	// Compose matrix and check result
	Matrix<double> U = L * Lt;// U = L*Lt

	U.Display();
*/


/*
	Matrix<double> A( 3, 3 );	// Matrix to be decomposed.
	Matrix<double> L( 3, 3 );	// Lower triangular matrix
	Matrix<double> D( 3, 3 );	// Lower triangular matrix

	A(0, 0) = 4;		A(0, 1) = 12;	A(0, 2) = -16;
	A(1, 0) = 12;	A(1, 1) = 37;	A(1, 2) = -43;
	A(2, 0) = -16;	A(2, 1) = -43;	A(2, 2) = 98;

	A.Display();

	//ModifiedCholesky( L, D, A );
	IncompleteCholesky( L, D, A );

	Matrix<double> Lt( 3, 3 );
	Transpose( Lt, L );

	L.Display();
	Lt.Display();

	auto validation = L * D * Lt;
	validation.Display();
*/


	// https://subscription.packtpub.com/book/big_data_and_business_intelligence/9781789346466/2/ch02lvl1sec21/the-cholesky-decomposition
	/*
	A =	[ 10, -1,  2,  0 ]
		[ -1, 11, -1,  3 ]
		[  2, -1, 10, -1 ]
		[  0,  3, -1,  8 ]

	x =	[ 1, 2, -1, 1 ]

	b =	[ 6, 25, -11, 15 ]
	*/

	DynamicMatrix<double> a(4, 4), b(4, 1), x(4, 1);

	a(0, 0) = 10;	a(0, 1) = -1;	a(0, 2) = 2;	a(0, 3) = 0;
	a(1, 0) = -1;	a(1, 1) = 11;	a(1, 2) = -1;	a(1, 3) = 3;
	a(2, 0) = 2;	a(2, 1) = -1;	a(2, 2) = 10;	a(2, 3) = -1;
	a(3, 0) = 0;	a(3, 1) = 3;	a(3, 2) = -1;	a(3, 3) = 8;

	b(0, 0) = 6;
	b(1, 0) = 25;
	b(2, 0) = -11;
	b(3, 0) = 15;

	tcout << "//======================== A ===========================//\n";
	a.Display();

	tcout << "//======================== b ===========================//\n";
	b.Display();


	tcout << "//==================== x(cholesky) =====================//\n";
	Cholesky_Solver<double>	cholesky_solver(a);
	Zero(x);
	cholesky_solver.Solve(x, b);
	x.Display();


	tcout << "//================ x(modified cholesky) ================//\n";
	ModifiedCholesky_Solver<double> mcholesky_solver(a);
	Zero(x);
	mcholesky_solver.Solve(x, b);
	x.Display();


	tcout << "//=============== x(incomplete cholesky) ===============//\n";
	IncompleteCholesky_Solver<double> icholesky_solver(a);
	Zero(x);
	icholesky_solver.Solve(x, b);
	x.Display();


	return 0;
}