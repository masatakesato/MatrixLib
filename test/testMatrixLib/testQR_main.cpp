#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/GramSchmidt.h>
#include	<matrixlib/QR.h>






int main()
{
	{
		tcout << _T("//============= QR_GramSchmidt ( 3x3 Matrix ) =============//\n");

		StaticMatrix<double, 3, 3> A;
		DynamicMatrix<double> /*v(3, 3),*/ Q(3, 3), R(3, 3);//, x(3, 3);

		A(0, 0) =0;
		A(1, 0) =1;
		A(2, 0) =1;

		A(0, 1) =1;
		A(1, 1) =0;
		A(2, 1) =1;

		A(0, 2) =1;
		A(1, 2) =1;
		A(2, 2) =0;

		

//		ModifiedGramSchmidt( v, q, x );
//		q.Display();
//		v.Display();

		QR_GramSchmidt( /*v,*/ Q, R, A );

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A - Q*R:\n");
		(A - Q*R).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//============= QR_GramSchmidt ( 4x3 Matrix ) =============//\n");

		DynamicMatrix<double> /*v(4, 3),*/ Q(4, 4), R(4, 3), A(4, 3);

		A(0, 0) = 1;
		A(1, 0) = 2;
		A(2, 0) = 1;
		A(3, 0) = 0;

		A(0, 1) = 1;
		A(1, 1) = 2;
		A(2, 1) = 1;//3;
		A(3, 1) = 0;//1;

		A(0, 2) = 1;
		A(1, 2) = 2;
		A(2, 2) = 3;
		A(3, 2) = 1;


		//ModifiedGramSchmidt( v, q, x );
		//q.Display();
		//v.Display();

		QR_GramSchmidt( /*v,*/ Q, R, A );

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A - Q*R:\n");
		(A - Q*R).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//============= QR_Householder ( 4x2 Matrix ) =============//\n");
		//Matrix<double> v(4, 3), Q(4, 4), R(4, 3), A(4, 3);

		//A(0, 0) = 1;
		//A(1, 0) = 1;
		//A(2, 0) = 1;
		//A(3, 0) = 1;

		//A(0, 1) = -1;
		//A(1, 1) = 4;
		//A(2, 1) = 4;
		//A(3, 1) = -1;

		//A(0, 2) = 4;
		//A(1, 2) = -2;
		//A(2, 2) = 2;
		//A(3, 2) = 0;

		DynamicMatrix<double> /*v(4, 1),*/ P(2, 2), Q(4, 4), R(4, 2), A(4, 2);
		//OreOreLib::Array<double>	rs(5);
		OreOreLib::Array<int>	Perm(2);

		A(0, 0) = 1;
		A(1, 0) = 2;
		A(2, 0) = 1;
		A(3, 0) = 0;

		A(0, 1) = 1;
		A(1, 1) = 2;
		A(2, 1) = 3;
		A(3, 1) = 1;

		QR_Householder( /*v,*/ Q, R, A );

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A - Q*R:\n");
		(A - Q*R).Display();
		tcout << tendl;



		tcout << _T("//============= QR_ColumnPivotHouseholder ( 4x2 Matrix ) =============//\n");
		//Zero(v);
		Zero(Q);
		Zero(R);

		int rank = QR_ColumnPivotHouseholder( /*v, rs,*/ Perm, Q, R, A );

		tcout << _T("Rank: ") << rank << tendl << tendl;

		tcout << _T("P:\n");
		TransposedPermutation( P, Perm );// column permutation
		P.Display();
		tcout << tendl;

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A*P - Q*R:\n");
		(A*P-Q*R).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{

		tcout << _T("//============= QR_ColumnPivotHouseholder ( 6x5 Matrix ) =============//\n");

		// TODO: http://faculty.nps.edu/borges/Teaching/MA3046/Matlab/qrlsdiary.html

		DynamicMatrix<double> /*v(6, 1),*/ P(5, 5), Q(6, 6)/*, R(6, 3), A(6, 3)*/;
		//StaticMatrix<double, 3, 3> P;
		StaticMatrix<double, 6, 5> R, A;
		//OreOreLib::Array<double>	rs(5);
		OreOreLib::Array<int>	Perm(5);


		A(0, 0) = 27;	A(0, 1) = 21;	A(0, 2) = 9;
		A(1, 0) = 25;	A(1, 1) = 27;	A(1, 2) = 18;
		A(2, 0) = 15;	A(2, 1) = 22;	A(2, 2) = 22;
		A(3, 0) = 2;	A(3, 1) = 7;	A(3, 2) = 29;
		A(4, 0) = 19;	A(4, 1) = 1;	A(4, 2) = 10;
		A(5, 0) = 12;	A(5, 1) = 22;	A(5, 2) = 7;
		
		A(0, 3) = 27;
		A(1, 3) = 25;
		A(2, 3) = 15;
		A(3, 3) = 2;
		A(4, 3) = 19;
		A(5, 3) = 12;

		A(0, 4) = 27;
		A(1, 4) = 25;
		A(2, 4) = 15;
		A(3, 4) = 2;
		A(4, 4) = 19;
		A(5, 4) = 12;
		

		int rank = QR_ColumnPivotHouseholder( /*v, rs,*/ Perm, Q, R, A );

		//for( int i=0; i<Perm.Length(); ++i )
		//	tcout << Perm[i] << ", ";
		//tcout << tendl;

		tcout << _T("Rank: ") << rank << tendl << tendl;

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("P:\n");
		TransposedPermutation(P, Perm);// column permutation
		P.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A*P - Q*R:\n");
		(A*P-Q*R).Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T("//==================Worse condition problem ===================//\n");

		//double a_[10][2] =
		//{
		//	{ 1, 5 },
		//	{ 1, 5 },
		//	{ 1, 6 },
		//	{ 1, 3 },
		//	{ 1, 5 },
		//	{ 1, 5 },
		//	{ 1, 5 },
		//	{ 1, 2 },
		//	{ 1, -5 },
		//	{ 1, 5 }
		//};

		double a_[10][2] =
		{
			//{ 5, 5 },
			//{ 5, 5 },
			//{ 6, 6 },
			//{ 3, 3 },
			//{ 5, 5 },
			//{ 5, 5 },
			//{ 5, 5 },
			//{ 2, 2 },
			//{ -5, -5 },
			//{ 5, 5 }

			{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 },
		{ 1, 1 }
		};


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


		DynamicMatrix<double>	P(2, 2), Q(10, 10), R(10, 2), A(10, 2, (double*)a_ );//x( 10, 2 ), b( 10, 2, (double *)b_ );
		OreOreLib::Array<int>	Perm(2);

		int rank = QR_ColumnPivotHouseholder( /*v, rs,*/ Perm, Q, R, A );

		//for( int i=0; i<Perm.Length(); ++i )
		//	tcout << Perm[i] << ", ";
		//tcout << tendl;

		tcout << _T("Rank: ") << rank << tendl << tendl;

		tcout << _T("A:\n");
		A.Display();
		tcout << tendl;

		tcout << _T("P:\n");
		TransposedPermutation(P, Perm);// column permutation
		P.Display();
		tcout << tendl;

		tcout << _T("Q:\n");
		Q.Display();
		tcout << tendl;

		tcout << _T("R:\n");
		R.Display();
		tcout << tendl;

		tcout << _T("A*P - Q*R:\n");
		(A*P-Q*R).Display();
		tcout << tendl;

	}

	return 0;
}