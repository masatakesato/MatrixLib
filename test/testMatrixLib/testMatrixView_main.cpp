#include <type_traits>
#include	<vector>
#include	<crtdbg.h>



#include	<matrixlib/MatrixLib.h>




//template< typename T, size_t... Args > class Mat; 
//
//
//
//
//
//// Dynamic Matrix
//template< typename T >
//class Mat<T> : Matrix<T>
//{
//public:
//
//	Mat()
//	{
//		tcout << "Dynamic Matrix\n";
//	}
//};
//
//
//// Static Matrix
//template< typename T, size_t Row, size_t Col >
//class Mat<T, Row, Col> : Matrix<T>
//{
//public:
//
//	Mat()
//	{
//		tcout << "Static Matrix\n";
//	}
//};
//
//
//// MatrixView
//template< typename T>
//class Mat<VIEW<T>> : Matrix<T>
//{
//public:
//
//	Mat()
//	{
//		tcout << "MatrixView\n";
//	}
//
//};
//
//Mat<double> matdynamic;
//Mat<double, 5, 5> matstatic;
//Mat<VIEW<double>> matview;




template< typename T >
void InitMatrix( IMatrix<T>& mat )
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

	Matrix<int, 30, 30>	sa;
	StaticMatrix<int, 30, 20>	sb;

	MatrixView<int> va(a, 0, 0, 3, 3);//MatrixView<int>	
	MatrixView<int> vb(b, 0, 0, 3, 2);//MatrixView<int>

	MatrixView<int> vsa(sa, 13, 26, 3, 3);//MatrixView<int>
	MatrixView<int> vsb(sb, 15, 15, 3, 2);//MatrixView<int>

	One(b);
	One(sb);

	{
		StaticMatrix<int, 3, 3> c(a);
		StaticMatrix<int, 3, 3> d(vsa);
		MatrixView<int> e(sa, 13, 26, 3, 3);
		e.Copy(vsa);
		e.Copy(c);
	}

	{
		std::vector< MatrixView<int> > vec;
		vec.push_back( MatrixView<int>(a, 0, 0, 3, 3) );
	}






	{
		tcout << "//============== Operators... ================//\n";

		InitMatrix(a);
		InitMatrix(vsa);
		InitMatrix(vsb);


		tcout << "vsa + a:\n";
		(vsa + a).Display();
		tcout << tendl;

		tcout << "a - vsa:\n";
		(a - vsa).Display();
		tcout << tendl;

		tcout << "vsa * a:\n";
		(vsa * a).Display();
		tcout << tendl;


		tcout << "a += vsa:\n";
		a += vsa;
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "vsa -= a:\n";
		vsa -= a;
		vsa.Display();
		tcout << tendl;
		InitMatrix(vsa);

		tcout << "vsa *= a;\n";
		vsa *= a;
		vsa.Display();
		tcout << tendl;
		InitMatrix(vsa);


		tcout << "a *= vsb;\n";
		auto tmp = a;
		tmp *= vsb;
		tmp.Display();
		tcout << tendl;
	}
	
	tcout << tendl;

	{
		tcout << "//============== Functions... ================//\n";

		InitMatrix(a);
		InitMatrix(vsa);
		InitMatrix(vsb);
		auto out = a;		


		tcout << "Add( out, vsa, a );\n";
		Add( out, vsa, a );
		out.Display();
		tcout << tendl;

		tcout << "Subtract( out, a, vsa );\n";
		Subtract( out, a, vsa );
		out.Display();
		tcout << tendl;

		tcout << "Multiply( out, vsa, a );\n";
		Multiply( out, vsa, a );
		out.Display();
		tcout << tendl;


		tcout << "Add( a, vsa );\n";
		Add( a, vsa );
		a.Display();
		tcout << tendl;
		InitMatrix(a);

		tcout << "Subtract( vsa, a );\n";
		Subtract( vsa, a );
		vsa.Display();
		tcout << tendl;
		InitMatrix(vsa);

		tcout << "Multiply( vsa, a );\n";
		Multiply( vsa, a );
		vsa.Display();
		tcout << tendl;
		InitMatrix(vsa);


		tcout << "Multiply( a, vsb );\n";
		auto tmp = a;
		Multiply( tmp, vsb );
		tmp.Display();
		tcout << tendl;
	}

	tcout << tendl;

	{
		tcout << "//============== ViewMatrix of ViewMatrix functions... ================//\n";
		One(sa);

		MatrixView<int> view(sa, 0, 0, sa.numRows(), sa.numCols());// vvsa_(vsa_, 1, 1, 26, 26);

		auto size = Min( sa.numRows(), sa.numCols() )-2;

		for( int i=size; i>=0; i-=2 )
		{
			if( i%4==0 )
				Zero(view);
			else
				One(view);
				
			view.Init( view, 1, 1, i, i );
		}

		sa.Display();

		tcout << tendl;
	}

	tcout << tendl;


	return 0;

}