#include	<matrixlib/GaussSeidel.h>


// http://pc-physics.com/jacobi1.html


int main( int argc, char *argv[] )
{
	tcout << "//===================== Test Gauss-Seidel method ===================//\n";

	/*
	A = [ 3, -6, 9 ]
		[ 2, 5, -8 ]
		[ 1, -4, 7 ]

	x = [ 3, 2, 1 ]

	b = [ 6, 8, 2 ]
	*/

	DynamicMatrix<float> A(3, 3), b(3, 1), x(3, 1);

	A(0, 0) = 3;		A(0, 1) = -6;	A(0, 2) = 9;
	A(1, 0) = 2;		A(1, 1) = 5;		A(1, 2) = -8;
	A(2, 0) = 1;		A(2, 1) = -4;	A(2, 2) = 7;

	b(0, 0) = 6;
	b(1, 0) = 8;
	b(2, 0) = 2;

	//x(0, 0) = 1;
	//x(1, 0) = 1;
	//x(2, 0) = 1;


	tcout << _T("A:\n");
	A.Display();
	tcout << tendl;

	tcout << _T("b:\n");
	b.Display();
	tcout << tendl;

	GaussSeidel( A, x, b, 50 );

	tcout << _T("x:\n");
	x.Display();
	tcout << tendl;

	tcout << _T("Ax - b:\n");
	(A*x - b).Display();
	tcout << tendl;

	return 0;
}