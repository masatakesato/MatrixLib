#include	<algorithm>

#include	<matrixlib/MatrixLib.h>



int main()
{
	DynamicMatrix<float> a(7, 7);
	StaticMatrix<float, 7, 7> sa;
	MatrixView<float> va(a, 0, 0, 7, 7);
	MatrixView<float> vsa(sa, 0, 0, 7, 7);	

	int perm[] = {0,1,2,3,4,5,6};
	int shuffled[] = {0,1,2,3,4,5,6};
	std::random_shuffle( &shuffled[0], &shuffled[6] );

	{
		tcout << _T( "//================ Matrix Permutation =================//\n" );

		Zero(a);
		a.Display();
		tcout << tendl;


		tcout << _T( "Init:\n" );
		Permutation( a, 7, perm );// row permutation
		a.Display();
		tcout << tendl;

		tcout << _T( "Shuffle:\n" );
		TransposedPermutation( a, 7, shuffled );// column permutation
		a.Display();
		tcout << tendl;

		tcout << _T( "Swap row 0 and 6:\n" );
		a.SwapRow(0, 6);
		a.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T( "//================ MatrixView Permutation =================//\n" );

		Zero(va);
		va.Display();
		tcout << tendl;


		tcout << _T( "Init:\n" );
		Permutation( va, 7, perm );// row permutation
		va.Display();
		tcout << tendl;

		tcout << _T( "Shuffle:\n" );
		TransposedPermutation( va, 7, shuffled );// column permutation
		va.Display();
		tcout << tendl;

		tcout << _T( "Swap row 0 and 6:\n" );
		va.SwapRow(0, 6);
		va.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T( "//================ StaticMatrix Permutation =================//\n" );

		Zero(sa);
		sa.Display();
		tcout << tendl;


		tcout << _T( "Init:\n" );
		Permutation( sa, 7, perm );// row permutation
		sa.Display();
		tcout << tendl;

		tcout << _T( "Shuffle:\n" );
		TransposedPermutation( sa, 7, shuffled );// column permutation
		a.Display();
		tcout << tendl;

		tcout << _T( "Swap row 0 and 6:\n" );
		sa.SwapRow(0, 6);
		sa.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << _T( "//================ StaticMatrixView Permutation =================//\n" );

		Zero(vsa);
		vsa.Display();
		tcout << tendl;


		tcout << _T( "Init:\n" );
		Permutation( vsa, 7, perm );// row permutation
		vsa.Display();
		tcout << tendl;

		tcout << _T( "Shuffle:\n" );
		TransposedPermutation( vsa, 7, shuffled );// column permutation
		vsa.Display();
		tcout << tendl;

		tcout << _T( "Swap row 0 and 6:\n" );
		vsa.SwapRow(0, 6);
		vsa.Display();
		tcout << tendl;
	}


	return 0;
}