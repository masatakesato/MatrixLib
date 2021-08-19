#include	<crtdbg.h>
#include	<vector>


#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/Substitution.h>



template< typename Type >
void InitMatrix( IMatrix<Type>& mat )
{
	for( int i=0; i<mat.numRows(); ++i )
		for( int j=0; j<mat.numCols(); ++j )
			mat(i, j) = i*mat.numCols() + j;
}




int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	DynamicMatrix<int> a(3, 3);
	DynamicMatrix<int> b(3, 2);
	StaticMatrix<int, 3, 3>	sa = { { -1, -2, -3 }, { 4, 5, 6 }, { -7, -8, -9 } };
	StaticMatrix<int, 3, 2>	sb;

	One(b);
	One(sb);

	//DynamicMatrix<float> c =
	//{
	//	{ -1, -2, -3.6f },
	//	{ 4,  5,  6  },
	//	{ -7, -8, -9 }
	//};

	//c.Display();
	//sa.Display();

	{
		tcout << "//============== Matrix operators... ================//\n";

		InitMatrix(a);

		tcout << "a:\n";
		a.Display();
		tcout << tendl;


		tcout << "a + a:\n";
		(a + a).Display();
		tcout << tendl;

		tcout << "a - a:\n";
		(a - a).Display();
		tcout << tendl;

		tcout << "a * a:\n";
		(a * a).Display();
		tcout << tendl;


		tcout << "a += a:\n";
		a += a;
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "a -= a:\n";
		a -= a;
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "a *= a:\n";
		a *= a;
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "a *= b:\n";
		auto c = a;
		c *= b;
		c.Display();
		tcout << tendl;

		tcout << "a *= sb:\n";
		c = a;
		c *= sb;
		c.Display();
		tcout << tendl;

	}

	tcout << tendl;

	{
		tcout << "//============== Matrix functions... ================//\n";

		auto out = a;
		InitMatrix(a);

		tcout << "a:\n";
		a.Display();
		tcout << tendl;


		tcout << "Add( out, a, a );\n";
		Add( out, a, a );
		out.Display();
		tcout << tendl;

		tcout << "Subtract( out, a, a );\n";
		Subtract( out, a, a );
		out.Display();
		tcout << tendl;

		tcout << "Multiply( out, a, a );\n";
		Multiply( out, a, a );
		out.Display();
		tcout << tendl;


		tcout << "Add( a, a );\n";
		Add( a, a );
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "Subtract( a, a );\n";
		Subtract( a, a );
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "Multiply( a, a );\n";
		Multiply( a, a );
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "Multiply( a, b );\n";
		auto c = a;
		Multiply( c, b );
		c.Display();
		tcout << tendl;

		tcout << "Multiply( a, sb );\n";
		c = a;
		Multiply( c, sb );
		c.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//============== StaticMatrix operators... ================//\n";

		InitMatrix(sa);

		tcout << "sa:\n";
		sa.Display();
		tcout << tendl;


		tcout << "sa + sa:\n";
		(sa + sa).Display();
		tcout << tendl;

		tcout << "sa - sa:\n";
		(sa - sa).Display();
		tcout << tendl;

		tcout << "sa * sa:\n";
		(sa * sa).Display();
		tcout << tendl;


		tcout << "sa += sa:\n";
		sa += sa;
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "sa -= sa:\n";
		sa -= sa;
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "sa *= sa:\n";
		sa *= sa;
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "sa *= a:\n";
		sa *= a;
		sa.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//============== StaticMatrix functions... ================//\n";

		auto out = sa;
		InitMatrix(sa);

		tcout << "sa:\n";
		sa.Display();
		tcout << tendl;


		tcout << "Add( out, sa, sa );\n";
		Add( out, sa, sa );
		out.Display();
		tcout << tendl;

		tcout << "Subtract( out, sa, sa );\n";
		Subtract( out, sa, sa );
		out.Display();
		tcout << tendl;

		tcout << "Multiply( out, sa, sa );\n";
		Multiply( out, sa, sa );
		out.Display();
		tcout << tendl;


		tcout << "Add( sa, sa );\n";
		Add( sa, sa );
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "Subtract( sa, sa );\n";
		Subtract( sa, sa );
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "Multiply( sa, sa );\n";
		Multiply( sa, sa );
		sa.Display();
		tcout << tendl;
		InitMatrix(sa);

		tcout << "Multiply( sa, a );\n";
		Multiply(sa, a);
		sa.Display();
		tcout << tendl;
	}

	tcout << tendl;


	return 0;
}