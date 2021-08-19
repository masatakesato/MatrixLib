#include	<matrixlib/GaussElimination.h>
#include	<matrixlib/LU.h>
#include	<matrixlib/Cholesky.h>






int main( int argc, char* argv )
{
	{
		tcout << _T("//=============== Inverse Matrix Calculation using Gaussian Elimination ================//\n" );
		DynamicMatrix<double> inv(3, 3), mat(3, 3);

		mat(0, 0) = 1;	mat(0, 1) = 2;	mat(0, 2) = 1;
		mat(1, 0) = 2;	mat(1, 1) = 1;	mat(1, 2) = 0;
		mat(2, 0) = 1;	mat(2, 1) = 1;	mat(2, 2) = 2;

		Inverse_Gauss( inv, mat );

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//=============== Inverse Matrix Calculation using LU decomposition ================//\n" );
		DynamicMatrix<double> inv(3, 3), mat(3, 3);

		mat(0, 0) = 1;	mat(0, 1) = 2;	mat(0, 2) = 1;
		mat(1, 0) = 2;	mat(1, 1) = 1;	mat(1, 2) = 0;
		mat(2, 0) = 1;	mat(2, 1) = 1;	mat(2, 2) = 2;

		Inverse_LU( inv, mat );

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//============= Inverse Matrix Calculation using LU_Pivot decomposition ============//\n" );
		DynamicMatrix<double> inv(3, 3), mat(3, 3);

		mat(0, 0) = 1;	mat(0, 1) = 2;	mat(0, 2) = 1;
		mat(1, 0) = 2;	mat(1, 1) = 1;	mat(1, 2) = 0;
		mat(2, 0) = 1;	mat(2, 1) = 1;	mat(2, 2) = 2;

		Inverse_LU_Pivot( inv, mat );

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//============= Inverse Matrix Calculation using LU_FullPivot decomposition ============//\n" );
		DynamicMatrix<double> inv(3, 3), mat(3, 3);

		mat(0, 0) = 1;	mat(0, 1) = 2;	mat(0, 2) = 1;
		mat(1, 0) = 2;	mat(1, 1) = 1;	mat(1, 2) = 0;
		mat(2, 0) = 1;	mat(2, 1) = 1;	mat(2, 2) = 2;

		Inverse_LU_FullPivot( inv, mat );

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		DynamicMatrix<double> mat(3, 3);

		mat(0, 0) = 1;	mat(0, 1) = 2;	mat(0, 2) = 1;
		mat(1, 0) = 2;	mat(1, 1) = 1;	mat(1, 2) = 0;
		mat(2, 0) = 1;	mat(2, 1) = 1;	mat(2, 2) = 2;


		tcout << _T("//=============== Inverse Matrix Calculation using LU_Solver ================//\n" );

		LU_Solver<double> solver(mat);
		auto inv = solver.Inverse();

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();	//solver.Inverse().Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;

		
		tcout << _T("//============ Inverse Matrix Calculation using PivotLU_Solver =============//\n" );

		PivotLU_Solver<double> psolver(mat);
		inv = psolver.Inverse();

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();	//solver.Inverse().Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
		
		
		tcout << _T("//============ Inverse Matrix Calculation using FullPivotLU_Solver =============//\n" );

		FullPivotLU_Solver<double> fpsolver(mat);
		fpsolver.SetCoeff(mat);
		inv = psolver.Inverse();

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();	//solver.Inverse().Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
		
	}

	tcout << tendl;

	{
		tcout << _T("//=============== Inverse Matrix Calculation using Cholesky decomposition ================//\n" );
		DynamicMatrix<double> inv(3, 3), mat(3, 3);

		mat(0, 0) = 4;		mat(0, 1) = 12;		mat(0, 2) = -16;
		mat(1, 0) = 12;		mat(1, 1) = 37;		mat(1, 2) = -43;
		mat(2, 0) = -16;	mat(2, 1) = -43;	mat(2, 2) = 98;

		Inverse_Cholesky( inv, mat );

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//=============== Inverse Matrix Calculation using Cholesky_Solver ================//\n" );
		DynamicMatrix<double> mat(3, 3);

		mat(0, 0) = 4;		mat(0, 1) = 12;		mat(0, 2) = -16;
		mat(1, 0) = 12;		mat(1, 1) = 37;		mat(1, 2) = -43;
		mat(2, 0) = -16;	mat(2, 1) = -43;	mat(2, 2) = 98;

		Cholesky_Solver<double>	solver(mat);
		auto& inv = solver.Inverse();

		tcout << _T("mat:\n");
		mat.Display();
		tcout << tendl;

		tcout << _T("mat^-1:\n");
		inv.Display();
		tcout << tendl;

		tcout << _T("mat * mat^-1:\n");
		(mat * inv).Display();
		tcout << tendl;
	}

	return 0;
}
