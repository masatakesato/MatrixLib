#include	<crtdbg.h>
#include	<vector>


#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/Substitution.h>




int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );



	DynamicMatrix<float> m;

	m.Init(3, 3);


	RowVector<float>	rowvec(2);
	rowvec.Init(66);


	ColVector<float>	colvec(4);
	colvec.Init(9);


	ColVector<float>	colvec2(6);

	colvec2 = rowvec;
	colvec2.Display();

	StaticRowVector<float, 4> srowvec;
	srowvec.Display();

	StaticColVector<float, 8>	scolvec;
	scolvec.Display();

	MatrixView< float> aa( m, 0, 0, 3, 3 );


	RowVectorView<float> v;
	v.Init( rowvec, 0, 0, 66);

	RowVectorView<float> vvv( rowvec, 0, 0, 2 );

	vvv(0)	= -9999.9f;
	vvv(1)	= 999999.0f;


//	vvv.Display();


//	rowvec.Display();


//	Matrix<VIEW<float>, DynamicDim, 2> staticview;

	//Matrix<VIEW<float>, 2, 2> staticmat;
	//MatrixView<float>



	// initializer_list test
	//RowVector<float> vec_row = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f };
	//ColVector<float> vec_col = { -1.1f, -2.2f, -3.3f, -4.4f, -5.5f };

	//vec_row.Display();
	//vec_col.Display();


	StaticRowVector<float, 5> svec_row = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f };
	StaticColVector<float, 5> svec_col = { -1.1f, -2.2f, -3.3f, -4.4f, -5.5f };

	svec_row.Display();
	svec_col.Display();


	return 0;
}