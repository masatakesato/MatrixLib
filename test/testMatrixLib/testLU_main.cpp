#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/LU.h>




int main( int argc, char *argv[] )
{	

	//###########################################################################################//

	{
		tcout << "//============================== LU =============================//\n";
	
		//double a_[4][4] = {
		//	{ 3, -7, -2, 2 },
		//	{ -3, 5, 1, 0 },
		//	{ 6, -4, 0, -5 },
		//	{ -9, 5, -5, 12 },
		//};
		
		//L = {
		//{  1,  0, 0, 0 },
		//{ -1,  1, 0, 0 },
		//{  2, -5, 1, 0 },
		//{ -3,  8, 3, 1 }
		//}

		//U = {
		//{ 3, -7, -2,  2 },
		//{ 0, -2, -1,  2 },
		//{ 0,  0, -1,  1 },
		//{ 0,  0,  0, -1 }
		//}

		double a_[]= {
			-1,	3,	0,	0,	4,	0,
			2,	-5,	0,	0,	0,	1,
			0,	0,	-2,	3,	0,	0,
			0,	0,	7,	-1,	0,	0,
			-3,	0,	0,	4,	6,	0,
			0,    5,    0,    0,   -7,    8
		};


		DynamicMatrix<double> L(6, 6), U(6, 6), A(6, 6, (double*)a_);

		tcout << "A:\n";
		A.Display();
		tcout << tendl;

		// Standard usage. 
		auto comp = A;
		LU( comp );

		// Alternative usage.
		//Matrix<double> comp(A.numRows(), A.numCols());
		//LU( comp, A );

		Zero(L);
		Zero(U);
		for( int i=0; i<comp.numRows(); ++i )
		{
			// L
			L(i, i) = 1;
			for( int j=0; j<i; ++j )
				L(i, j) = comp(i, j);

			// U
			for( int j=i; j<comp.numCols(); ++j )
				U(i, j) = comp(i, j);
		}

		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		tcout << "A-L*U:\n";
		(A-L*U).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//=========================== LU_Pivot ==========================//\n";

		double a_[4][4] = {
			{ 3, -7, -2, 2 },
			{ -3, 5, 1, 0 },
			{ 6, -4, 0, -5 },
			{ -9, 5, -5, 12 },
		};

		//L = {
		//{  1,  0, 0, 0 },
		//{ -1,  1, 0, 0 },
		//{  2, -5, 1, 0 },
		//{ -3,  8, 3, 1 }
		//}

		//U = {
		//{ 3, -7, -2,  2 },
		//{ 0, -2, -1,  2 },
		//{ 0,  0, -1,  1 },
		//{ 0,  0,  0, -1 }
		//}
		//


		DynamicMatrix<double> L(4, 4), U(4, 4), A(4, 4, (double*)a_);

		OreOreLib::Array<int> pivot( Min(A.numRows(), A.numCols()) );
	
		tcout << "A:\n";
		A.Display();
		tcout << tendl;

		// Standard usage. 
		auto comp = A;
		LU_Pivot( comp, pivot );

		// Alternative usage.
		//Matrix<double> comp(A.numRows(), A.numCols());
		//LU_Pivot( comp, pivot, A );

		Zero(L);
		Zero(U);
		for( int i=0; i<comp.numRows(); ++i )
		{
			// L
			L(i, i) = 1;
			for( int j=0; j<i; ++j )
				L(i, j) = comp(i, j);

			// U
			for( int j=i; j<comp.numCols(); ++j )
				U(i, j) = comp(i, j);
		}

		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		DynamicMatrix<double> P(4, 4);
		Permutation( P, pivot );// row permutation

		tcout << "P*A - L*U:\n";
		(P*A - L*U).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//=========================== LU_FullPivot ==========================//\n";

		double a_[4][4] = {
			{ 3, -7, -2, 2 },
			{ -3, 5, 1, 0 },
			{ 6, -4, 0, -5 },
			{ -9, 5, -5, 12 },
		};

		
		//L = {
		//{  1,  0, 0, 0 },
		//{ -1,  1, 0, 0 },
		//{  2, -5, 1, 0 },
		//{ -3,  8, 3, 1 }
		//}

		//U = {
		//{ 3, -7, -2,  2 },
		//{ 0, -2, -1,  2 },
		//{ 0,  0, -1,  1 },
		//{ 0,  0,  0, -1 }
		//}
		//

		DynamicMatrix<double> L(4, 4), U(4, 4), A(4, 4, (double*)a_);

		OreOreLib::Array<int> p( A.numRows() ), q( A.numCols() );
	
		tcout << "A:\n";
		A.Display();
		tcout << tendl;

		// Standard usage. 
		auto comp = A;
		LU_FullPivot( comp, p, q );

		// Alternative usage.
		//Matrix<double> comp(A.numRows(), A.numCols());
		//LU_FullPivot( comp, pivot, A );

		Zero(L);
		Zero(U);
		for( int i=0; i<comp.numRows(); ++i )
		{
			// L
			L(i, i) = 1;
			for( int j=0; j<i; ++j )
				L(i, j) = comp(i, j);

			// U
			for( int j=i; j<comp.numCols(); ++j )
				U(i, j) = comp(i, j);
		}

		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		DynamicMatrix<double> P(4, 4);
		TransposedPermutation( P, p );// Inverse of row permutation

		DynamicMatrix<double> Q(4, 4);
		Permutation( Q, q );// Inverse of col permutation = normal permutation

		tcout << "P:\n";
		p.Display();
		tcout << tendl;

		tcout << "Q:\n";
		q.Display();
		tcout << tendl;

		tcout << "L*U:\n";
		(L*U).Display();
		tcout << tendl;

		tcout << "A - P*L*U*Q:\n";
		(A - P*L*U*Q).Display();
		tcout << tendl;
	}
	
	tcout << tendl;

	{
		tcout << "//=========================== ILU(0) ==========================//\n";
	

		// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.330.2112&rep=rep1&type=pdf 

		double a_[4][4] = {
			{ 1, 0, 3, 4 },
			{ 2, 1, 0, 6 },
			{ 0, 0, 2, 1 },
			{ 3, 0, 5, 3 },
		};

		// LU result:
		//L = {
		//	{ 1,  0,  0,  0 },
		//	{ 2,  1,  0,  0 },
		//	{ 0,  0,  1,  0 },
		//	{ 3,  0, -2,  1 }
		//}
		//U = {
		//	{ 1,  0,  3,  4 },
		//	{ 0,  1, -6, -2 },
		//	{ 0,  0,  2,  1 },
		//	{ 0,  0,  0, -7 }
		//}

		// ILU(0) result:
		//L = {
		//	{ 1,  0,  0,  0 },
		//	{ 2,  1,  0,  0 },
		//	{ 0,  0,  1,  0 },
		//	{ 3,  0, -2,  1 }
		//}
		//U = {
		//	{ 1,  0,  3,  4 },
		//	{ 0,  1,  0, -2 },
		//	{ 0,  0,  2,  1 },
		//	{ 0,  0,  0, -7 }
		//}

		DynamicMatrix<double> L(4, 4), U(4, 4), A(4, 4, (double*)a_);

		tcout << "A:\n";
		A.Display();
		tcout << tendl;

		// Standard usage. 
		auto comp = A;
		ILU0( comp );

		// Alternative usage.
		//Matrix<double> comp(A.numRows(), A.numCols());
		//ILU0( comp, A );

		Zero(L);
		Zero(U);
		for( int i=0; i<comp.numRows(); ++i )
		{
			// L
			L(i, i) = 1;
			for( int j=0; j<i; ++j )
				L(i, j) = comp(i, j);

			// U
			for( int j=i; j<comp.numCols(); ++j )
				U(i, j) = comp(i, j);
		}

		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		tcout << "R(A-L*U):\n";
		(A-L*U).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//=========================== MILU(0) ==========================//\n";


		// http://forge.scilab.org/index.php/p/spilu/source/tree/HEAD/tests/unit_tests/milu0.tst

		double a_[]= {
			-1,	3,	0,	0,	4,	0,
			2,	-5,	0,	0,	0,	1,
			0,	0,	-2,	3,	0,	0,
			0,	0,	7,	-1,	0,	0,
			-3,	0,	0,	4,	6,	0,
			0,    5,    0,    0,   -7,    8
		};

		//expected L = {
		//	1,		0,			0,			0,			0,			0,
		//	-2,		1,			0,			0,			0,			0,
		//	0,		0,			1,			0,			0,			0,
		//	0,		0,			-3.5,		1,			0,			0,
		//	3,		0,			0,			0.4210526,	1,			0,
		//	0,		0.5555556,	0,			0,			0.4666667,	1
		// };

		//expected U = {
		//	-1,		3,			0,			0,			4,			0,
		//	0,		9,			0,			0,			0,			1,
		//	0,		0,			-2,			3,			0,			0,
		//	0,		0,			0,			9.5,		0,			0,
		//	0,		0,			0,			0,			-15,		0,
		//	0,		0,			0,			0,			0,			7.4444444
		//};


		DynamicMatrix<double> L(6, 6), U(6, 6), A(6, 6, (double*)a_);

		tcout << "A:\n";
		A.Display();
		tcout << tendl;

		// Standard usage. 
		auto comp = A;
		MILU0( comp );

		// Alternative usage.
		//Matrix<double> comp(A.numRows(), A.numCols());
		//MILU0( comp, A );

		Zero(L);
		Zero(U);
		for( int i=0; i<comp.numRows(); ++i )
		{
			// L
			L(i, i) = 1;
			for( int j=0; j<i; ++j )
				L(i, j) = comp(i, j);

			// U
			for( int j=i; j<comp.numCols(); ++j )
				U(i, j) = comp(i, j);
		}

		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		tcout << "R(=A-L*U):\n";
		(A-L*U).Display();
		tcout << tendl;
	}

	tcout << tendl; 

	{
		tcout << "//=========================== ILUT ==========================//\n";


		double a_[4][4] = {
			{ 3, -7, -2, 2 },
			{ -3, 5, 1, 0 },
			{ 6, -4, 0, -5 },
			{ -9, 5, -5, 12 },
		};
		
		//L = {
		//{  1,  0, 0, 0 },
		//{ -1,  1, 0, 0 },
		//{  2, -5, 1, 0 },
		//{ -3,  8, 3, 1 }
		//}

		//U = {
		//{ 3, -7, -2,  2 },
		//{ 0, -2, -1,  2 },
		//{ 0,  0, -1,  1 },
		//{ 0,  0,  0, -1 }
		//}
		//

		DynamicMatrix<double> L(4, 4), U(4, 4), A(4, 4, (double*)a_);

		tcout << "A:\n";
		A.Display();
		tcout << tendl;


		Zero(L);
		Zero(U);

		// Standard usage.
		auto comp = A;
		ILUT( L, U, comp, /*1.0e-2*/0.0 );


		tcout << "L:\n";
		L.Display();
		tcout << tendl;

		tcout << "U:\n";
		U.Display();
		tcout << tendl;

		tcout << "R(=A-L*U):\n";
		(A-L*U).Display();
		tcout << tendl;

	}
		
	return 0;
}