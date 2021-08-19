#include	<matrixlib/Solver.h>
#include	<matrixlib/LU.h>
#include	<matrixlib/Cholesky.h>
#include	<matrixlib/QR.h>



int main( int argc, char *argv[] )
{
	{
		StaticMatrix<double, 4, 4> A_;
		MatrixView<double>	A(A_, 0, 0, 4, 4);
		DynamicMatrix<double> b(4, 1), x(4, 1);

		A(0, 0) = 8;	A(0, 1) = 16;	A(0, 2) = 24;	A(0, 3) = 32;
		A(1, 0) = 2;	A(1, 1) = 7;	A(1, 2) = 12;	A(1, 3) = 17;
		A(2, 0) = 6;	A(2, 1) = 17;	A(2, 2) = 32;	A(2, 3) = 59;
		A(3, 0) = 7;	A(3, 1) = 22;	A(3, 2) = 46;	A(3, 3) = 105;

		b(0, 0) = 160;
		b(1, 0) = 70;
		b(2, 0) = 198;
		b(3, 0) = 291;

		//x(0, 0) = 4;
		//x(1, 0) = 3;
		//x(2, 0) = 2;
		//x(3, 0) = 1;

		//b(0, 0) = 0;
		//b(1, 0) = 0;
		//b(2, 0) = 0;
		//b(3, 0) = 0;


		LU_Solver<double> lu_solver(A);
		PivotLU_Solver<double>	pivotlu_solver(A);
		FullPivotLU_Solver<double>	fullpivotlu_solver(A);
		ILU0_Preconditioner<double>	ilu0_solver(A);
		MILU0_Preconditioner<double> milu0_solver(A);
		QRGramSchmidt_Solver<double> qrgs_solver(A);
		QRHouseholder_Solver<double> qrh_solver(A);
		QRColumnPivotHouseholder_Solver<double> qrcp_solver(A);

		tcout << "//==================== LU_Solver ==================//\n";
		Zero(x);
		lu_solver.Solve(x, b);
		x.Display();

		tcout << "//================== PivotLU_Solver ================//\n";
		Zero(x);
		pivotlu_solver.Solve(x, b);
		x.Display();

		tcout << "//================ FullPivotLU_Solver ==============//\n";
		Zero(x);
		fullpivotlu_solver.Solve(x, b);
		x.Display();
		
		tcout << "//================ ILU0_Preconditioner ==============//\n";
		Zero(x);
		ilu0_solver.Solve(x, b);
		x.Display();
		
		tcout << "//================ MILU0_Preconditioner ==============//\n";
		Zero(x);
		milu0_solver.Solve(x, b);
		x.Display();

		tcout << "//=============== QR GramSchmidt Solver ==============//\n";
		Zero(x);
		qrgs_solver.Solve(x, b);
		x.Display();

		tcout << "//=============== QR Householder Solver ==============//\n";
		Zero(x);
		qrh_solver.Solve(x, b);
		x.Display();

		tcout << "//======== QR Column Pivot Householder Solver ========//\n";
		Zero(x);
		qrcp_solver.Solve(x, b);
		tcout << "rank = " << qrcp_solver.Rank() << tendl;
		x.Display();

	}

	tcout << tendl;

	{
		DynamicMatrix<double> A(4, 4), b(4, 1), x(4, 1);

		A(0, 0) = 10;	A(0, 1) = -1;	A(0, 2) = 2;	A(0, 3) = 0;
		A(1, 0) = -1;	A(1, 1) = 11;	A(1, 2) = -1;	A(1, 3) = 3;
		A(2, 0) = 2;	A(2, 1) = -1;	A(2, 2) = 10;	A(2, 3) = -1;
		A(3, 0) = 0;	A(3, 1) = 3;	A(3, 2) = -1;	A(3, 3) = 8;

		b(0, 0) = 6;
		b(1, 0) = 25;
		b(2, 0) = -11;
		b(3, 0) = 15;

		//x(0, 0) = 1;
		//x(1, 0) = 2;
		//x(2, 0) = -1;
		//x(3, 0) = 1;


		Cholesky_Solver<double>	cholesky_solver(A);
		ModifiedCholesky_Solver<double> mcholesky_solver(A);
		IncompleteCholesky_Solver<double> icholesky_solver(A);
		
		tcout << "//=================== Cholesky_Solver =================//\n";
		Zero(x);
		cholesky_solver.Solve(x, b);
		x.Display();

		tcout << "//============== Modified Cholesky_Solver =============//\n";
		Zero(x);
		mcholesky_solver.Solve(x, b);
		x.Display();

		tcout << "//============= Incomplete Cholesky_Solver =============//\n";
		Zero(x);
		icholesky_solver.Solve(x, b);
		x.Display();
	}

	tcout << tendl;

	{
		DynamicMatrix<double> A(4, 4), b(4, 3), x(4, 3);

		A(0, 0) = 8;	A(0, 1) = 16;	A(0, 2) = 24;	A(0, 3) = 32;
		A(1, 0) = 2;	A(1, 1) = 7;	A(1, 2) = 12;	A(1, 3) = 17;
		A(2, 0) = 6;	A(2, 1) = 17;	A(2, 2) = 32;	A(2, 3) = 59;
		A(3, 0) = 7;	A(3, 1) = 22;	A(3, 2) = 46;	A(3, 3) = 105;

		b(0, 0) = 160;	b(0, 1) = -160;	b(0, 2) = 144;
		b(1, 0) = 70;	b(1, 1) = -70;	b(1, 2) = 84;
		b(2, 0) = 198;	b(2, 1) = -198;	b(2, 2) = 320;
		b(3, 0) = 291;	b(3, 1) = -291;	b(3, 2) = 615;

		//x(0, 0) = 4;
		//x(1, 0) = 3;
		//x(2, 0) = 2;
		//x(3, 0) = 1;

		//x(0, 1) = -4;
		//x(1, 1) = -3;
		//x(2, 1) = -2;
		//x(3, 1) = -1;

		//x(0, 2) = -5;
		//x(1, 2) = 6;
		//x(2, 2) = -7;
		//x(3, 2) = 8;


		LU_Solver<double> lu_solver(A, 2);
		PivotLU_Solver<double> pivotlu_solver(A, 3);
		FullPivotLU_Solver<double> fullpivotlu_solver(A, 3);
		ILU0_Preconditioner<double>	ilu0_solver(A, 3);
		MILU0_Preconditioner<double> milu0_solver(A, 3);
		ILUT_Preconditioner<double> ilut_precond(A, 0.0, 3);
		QRGramSchmidt_Solver<double> qrgs_solver(A);
		QRHouseholder_Solver<double> qrh_solver(A);
		QRColumnPivotHouseholder_Solver<double> qrcp_solver(A);

		tcout << "//==================== LU_Solver with multiple RHS ==================//\n";
		Zero(x);
		lu_solver.Solve(x, b);
		//lu_solver.Solve_Parallel(x, b);// Available if openmp is enabled.
		x.Display();

		tcout << "//================== PivotLU_Solver with multiple RHS ================//\n";
		Zero(x);
		pivotlu_solver.Solve(x, b);
		//pivotlu_solver.Solve_Parallel(x, b);// Available if openmp is enabled.
		x.Display();

		tcout << "//================ FullPivotLU_Solver with multiple RHS ==============//\n";
		Zero(x);
		fullpivotlu_solver.Solve(x, b);
		//pivotlu_solver.Solve_Parallel(x, b);// Available if openmp is enabled.
		x.Display();

		tcout << "//================ ILU0_Preconditioner with multiple RHS ==============//\n";
		Zero(x);
		ilu0_solver.Solve(x, b);
		x.Display();

		tcout << "//================ MILU0_Preconditioner with multiple RHS ==============//\n";
		Zero(x);
		milu0_solver.Solve(x, b);
		x.Display();

		tcout << "//================ ILUT_Preconditioner with multiple RHS ==============//\n";
		Zero(x);
		ilut_precond.Solve(x, b);
		x.Display();

		tcout << "//=============== QR GramSchmidt Solver with multiple RHS ==============//\n";
		Zero(x);
		qrgs_solver.Solve(x, b);
		x.Display();

		tcout << "//=============== QR Householder Solver with multiple RHS ==============//\n";
		Zero(x);
		qrh_solver.Solve(x, b);
		x.Display();

		tcout << "//======== QR Column Pivot Householder Solver with multiple RHS ========//\n";
		Zero(x);
		qrcp_solver.Solve(x, b);
		tcout << "rank = " << qrcp_solver.Rank() << tendl;
		x.Display();


	}

	tcout << tendl;

	{
		DynamicMatrix<double> A(4, 4), b(4, 3), x(4, 3);

		A(0, 0) = 10;	A(0, 1) = -1;	A(0, 2) = 2;	A(0, 3) = 0;
		A(1, 0) = -1;	A(1, 1) = 11;	A(1, 2) = -1;	A(1, 3) = 3;
		A(2, 0) = 2;	A(2, 1) = -1;	A(2, 2) = 10;	A(2, 3) = -1;
		A(3, 0) = 0;	A(3, 1) = 3;	A(3, 2) = -1;	A(3, 3) = 8;

		b(0, 0) = 6;
		b(1, 0) = 25;
		b(2, 0) = -11;
		b(3, 0) = 15;

		b(0, 1) = 33;
		b(1, 1) = 28;
		b(2, 1) = -14;
		b(3, 1) = 3;

		b(0, 2) = 37;
		b(1, 2) = -129;
		b(2, 2) = 61;
		b(3, 2) = -95;

		//x(0, 0) = 1;
		//x(1, 0) = 2;
		//x(2, 0) = -1;
		//x(3, 0) = 1;

		//x(0, 1) = 4;
		//x(1, 1) = 3;
		//x(2, 1) = -2;
		//x(3, 1) = -1;

		//x(0, 2) = 2;
		//x(1, 2) = -9;
		//x(2, 2) = 4;
		//x(3, 2) = -8;


		Cholesky_Solver<double>	cholesky_solver(A, 2);
		ModifiedCholesky_Solver<double> mcholesky_solver(A, 3);
		IncompleteCholesky_Solver<double> icholesky_solver(A, 4);

		tcout << "//=================== Cholesky_Solver multiple RHS =================//\n";
		Zero(x);
		cholesky_solver.Solve(x, b);
		x.Display();

		tcout << "//============== Modified Cholesky_Solver multiple RHS =============//\n";
		Zero(x);
		mcholesky_solver.Solve(x, b);
		x.Display();

		tcout << "//============= Incomplete Cholesky_Solver multiple RHS =============//\n";
		Zero(x);
		icholesky_solver.Solve(x, b);
		x.Display();
	}


	return 0;
}