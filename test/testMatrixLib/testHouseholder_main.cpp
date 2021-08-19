#include	<matrixlib/Householder.h>



int main()
{
	tcout << _T("//=============== Householder transformation ================//\n" );

	DynamicMatrix<double> A(4, 4), v(4, 1);

	v(0, 0) = -1;
	v(1, 0) = 1;
	v(2, 0) = 1;
	v(3, 0) = 1;


	tcout <<_T("v:\n");
	v.Display();
	tcout <<_T("\n");

	Householder( A, v );

	tcout <<_T("Transform result:\n");
	A.Display();

	return 0;
}
