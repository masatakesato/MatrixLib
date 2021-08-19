#include	<matrixlib/LU.h>
#include	<matrixlib/CG.h>
#include	<matrixlib/PCG.h>
#include	<matrixlib/Cholesky.h>


// http://pc-physics.com/jacobi1.html


void func( double *a )
{

}


int main( int argc, char *argv[] )
{
	tcout << "Test Conjugate Gradient method\n";

	
	//double a_[3][3] = {
	//	{  1,  2, -3 },
	//	{  2,  5, -4 },
	//	{ -3, -4,  8 } };

	////double x_[] = { 0, 1, 2 };

	//double b_[] = { -4, -3, 12 };
	//

	//Matrix<double>	A( 3, 3, (double *)a_ ),
	//				x( 3, 1 ),
	//				b( 3, 1, (double *)b_ );


	//double a_[10][10] = {
	//	{ 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	//	{ 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	//	{ 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	//	{ 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	//	{ 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0 },
	//	{ 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0 },
	//	{ 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0 },
	//	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0 },
	//	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0 },
	//	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0 } };

	//double x_[] = { 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 };

	//double b_[] = { 3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0 };


	// Good condition problem
	//double a_[10][10] =
	//{
	//	{ 94,  0,   0,   0,    0,   28,  0,   0,   32,  0 },
	//	{ 0,   59,  13,  5,    0,   0,   0,   10,  0,   0 },
	//	{ 0,   13,  72,  34,   2,   0,   0,   0,   0,   65 },
	//	{ 0,   5,   34,  114,  0,   0,   0,   0,   0,   55 },
	//	{ 0,   0,   2,   0,    70,  0,   28,  32,  12,  0 },
	//	{ 28,  0,   0,   0,    0,   87,  20,  0,   33,  0 },
	//	{ 0,   0,   0,   0,    28,  20,  71,  39,  0,   0 },
	//	{ 0,   10,  0,   0,    32,  0,   39,  46,  8,   0 },
	//	{ 32,  0,   0,   0,    12,  33,  0,   8,   82,  11 },
	//	{ 0,   0,   65,  55,   0,   0,   0,   0,   11,  100 } };

	double b_[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };//

	// Worse condition problem
	double a_[10][10] =
	{
		{ 894,    0,   0,     0,   0,   28,  0,   0,   1000,  70000 },
		{ 0,      5,   13,    5,   0,   0,   0,   0,   0,     0 },
		{ 0,      13,  72,    34,  0,   0,   0,   0,   0,     6500 },
		{ 0,      5,   34,    1,   0,   0,   0,   0,   0,     55 },
		{ 0,      0,   0,     0,   70,  0,   28,  32,  12,    0 },
		{ 28,     0,   0,     0,   0,   87,  20,  0,   33,    0 },
		{ 0,      0,   0,     0,   28,  20,  71,  39,  0,     0 },
		{ 0,      0,   0,     0,   32,  0,   39,  46,  8,     0 },
		{ 1000,   0,   0,     0,   12,  33,  0,   8,   82,    11 },
		{ 70000,  0,   6500,  55,  0,   0,   0,   0,   11,    100 } };


	//0.0174039
	//1.55702
	//-0.179325
	//-0.890778
	//0.0458074
	//0.13817
	//-0.119308
	//0.139099
	//-0.276422
	//0.00368562


	DynamicMatrix<double>	A( 10, 10, (double *)a_ );
	DynamicMatrix<double>	x( 10, 1 );//StaticMatrix<double, 10, 1>	x;//
	DynamicMatrix<double>	b( 10, 1, (double *)b_);//StaticMatrix<double, 10, 1>	b( b_);//


	// Preconditioner
	LU_Solver<double>					lu_precond(A);
	PivotLU_Solver<double>				pivotlu_precond(A);
	ILU0_Preconditioner<double>			ilu0_precond(A);
	Cholesky_Solver<double>				cholesky_precond(A);
	IncompleteCholesky_Solver<double>	icholesky_precond(A);


	tcout << "//============================== CG =============================//\n";
	One(x);//Random(x, -999.9, 999.9 );
	CG( A, x, b );
	x.Display();
	
	/*
	// Deprecated
	tcout << "//============================= ICCG ============================//\n";
	//Random(x, -999.9, 999.9 );
	One(x);//Zero(x);
	ICCG( A, x, b );
	x.Display();
	*/


	tcout << "//=========================== PCG(LU) ===========================//\n";
	One(x);
	PCG( A, x, b, lu_precond );
	x.Display();


	tcout << "//=========================== PCG(LUPivot) ===========================//\n";
	One(x);
	PCG( A, x, b, pivotlu_precond );
	x.Display();


	tcout << "//=========================== PCG(ILU(0)) ===========================//\n";
	One(x);
	PCG( A, x, b, ilu0_precond );
	x.Display();
	

	tcout << "//======================== PCG(Cholesky) ========================//\n";
	One(x);
	PCG( A, x, b, cholesky_precond );
	x.Display();


	tcout << "//=================== PCG(Incomplete Cholesky) ===================//\n";
	One(x);
	PCG( A, x, b, icholesky_precond );
	x.Display();



	tcout << "//========================== BiCG ===============================//\n";
	One(x);//Random(x, -999.9, 999.9 );
	BiCG( A, x, b );
	x.Display();


	/*
	// Deprecated.
	tcout << "//========================= ICBiCG ==============================//\n";
	//Random(x, -999.9, 999.9 );
	One(x);//Random(x);
	ICBiCG( A, x, b );
	x.Display();
	*/


	tcout << "//======================== PBiCG(LU) ===========================//\n";
	One(x);//Random(x);
	PBiCG( A, x, b, lu_precond );
	x.Display();


	tcout << "//======================== PBiCG(PivotLU) ===========================//\n";
	One(x);//Random(x);
	PBiCG( A, x, b, pivotlu_precond );
	x.Display();


	tcout << "//====================== PBiCG(ILU(0)) =========================//\n";
	One(x);//Random(x);
	PBiCG( A, x, b, ilu0_precond );
	x.Display();


	tcout << "//===================== PBiCG(Cholesky) ========================//\n";
	One(x);//Random(x);
	PBiCG( A, x, b, cholesky_precond );
	x.Display();


	tcout << "//================== PBiCG(Incomplete Cholesky) ================//\n";
	One(x);//Random(x);
	PBiCG( A, x, b, icholesky_precond );
	x.Display();



	tcout << "//========================= BiCGSTAB ==============================//\n";
	One(x);//Random(x, -999.9, 999.9 );
	BiCGSTAB( A, x, b );
	x.Display();


	/*
	// Deprecated.
	tcout << "//========================= ICBiCGSTAB ===========================//\n";
	//Random(x, -999.9, 999.9 );
	//Random( x );
	Zero(x);//One(x);//
	ICBiCGSTAB( A, x, b );
	x.Display();
	*/


	tcout << "//======================== PBiCGSTAB(LU) =========================//\n";
	Zero(x);//One(x);//Random( x );
	PBiCGSTAB( A, x, b, lu_precond );
	x.Display();


	tcout << "//======================== PBiCGSTAB(PivotLU) =========================//\n";
	Zero(x);//One(x);//Random( x );
	PBiCGSTAB( A, x, b, pivotlu_precond );
	x.Display();


	tcout << "//====================== PBiCGSTAB(ILU(0)) ========================//\n";
	One(x);//Random(x, -999.9, 999.9 );
	PBiCGSTAB( A, x, b, ilu0_precond );
	x.Display();


	tcout << "//====================== PBiCGSTAB(Cholesky) ========================//\n";
	One(x);//Random(x, -999.9, 999.9 );
	PBiCGSTAB( A, x, b, cholesky_precond );
	x.Display();


	tcout << "//================ PBiCGSTAB(Incomplete Cholesky) ================//\n";
	Zero(x);//One(x);//Random( x );
	PBiCGSTAB( A, x, b, icholesky_precond );
	x.Display();



	return 0;
}


/*

// http://www-section.cocolog-nifty.com/blog/2008/12/bicgstab-30b8-1.html

Matrix a = { {5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0} };

Vector b = { 3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0 };


Vector x = { 1, -1, 2, -2, 3, -3, 4, -4, 5, -5 };

*/

