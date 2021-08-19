#include	<matrixlib/MatrixLib.h>



int main()
{

	{
		tcout << _T( "//================= 3x2 Matrix Transpose test =================//\n" );

		DynamicMatrix<double> buf(20, 20);
		DynamicMatrix<double> /*c(3, 2),*/ ct(2, 3)/*, out22(2, 2), out33(3, 3)*/;
		StaticMatrix<double, 3, 2> c; 
		
		MatrixView<double> out33(buf, 1, 2, 3, 3), out22(buf, 14, 17, 2, 2);

		c(0,0) = 1;		c(0,1) = -3;
		c(1,0) = 5;		c(1,1) = -7;
		c(2,0) = -9;	c(2,1) = 11;

		Transpose(ct, c);

		tcout << _T( "c:\n" );
		c.Display();
		tcout << tendl;

		tcout << _T( "c^t:\n" );
		ct.Display();
		tcout << tendl;

		//  c^t * c =
		//   [  107   137 ]
		//   [ -137   179 ]

		tcout << _T( "TransposeMultiply(out, c, c);\n" );
		Zero(out22);
		TransposeMultiply(out22, c, c);
		out22.Display();
		tcout << tendl;

		tcout << _T( "TransposeMultiply(out, c);\n" );
		Zero(out22);
		TransposeMultiply(out22, c);
		out22.Display();
		tcout << tendl;

		tcout << _T( "Transpose(ct, c) and Multiply(out, ct, c);\n" );
		Zero(out22);
		Multiply(out22, ct, c);
		out22.Display();
		tcout << tendl;


		//  c * c^t =
		//   [  10   26   -42  ]
		//   [  26   74   -122 ]
		//   [ -42  -122   202 ]

		tcout << _T( "MultiplyTranspose(out, c, c);\n" );
		Zero(out33);
		MultiplyTranspose(out33, c, c);
		out33.Display();
		tcout << tendl;

		tcout << _T( "MultiplyTranspose(out, c);\n" );
		Zero(out33);
		MultiplyTranspose(out33, c);
		out33.Display();
		tcout << tendl;

		tcout << _T( "Transpose(ct, c) and Multiply(out, c, ct);\n" );
		Zero(out33);
		Multiply(out33, c, ct);
		out33.Display();
		tcout << tendl;

		buf.Display();
	}

	tcout << tendl;

	{
		tcout << _T( "//================= 2x3 Matrix Transpose test =================//\n" );

		DynamicMatrix<double> c(2, 3);//, /*ct(3, 2),*/ out33(3, 3), out22(2, 2);
		StaticMatrix<double, 3, 2> ct; 
		DynamicMatrix<double> buf(20, 20);
		MatrixView<double> out33(buf, 1, 2, 3, 3), out22(buf, 14, 17, 2, 2);


		c(0,0) = 1;		c(0,1) = 5;		c(0,2) = -9;	
		c(1,0) = -3;	c(1,1) = -7;	c(1,2) = 11;

		Transpose(ct, c);

		tcout << _T( "c:\n" );
		c.Display();
		tcout << tendl;

		tcout << _T( "c^t:\n" );
		ct.Display();
		tcout << tendl;


		//  c^t * c =
		//   [  10   26   -42  ]
		//   [  26   74   -122 ]
		//   [ -42  -122   202 ]

		tcout << _T( " TransposeMultiply(out, c, c);\n" );
		Zero(out33);
		TransposeMultiply(out33, c, c);
		out33.Display();

		tcout << _T( "TransposeMultiply(out,  c);\n" );
		Zero(out33);
		TransposeMultiply(out33, c);
		out33.Display();

		tcout << _T( "Transpose(ct, c) and Multiply(out, ct, c);\n" );
		Zero(out33);
		Transpose(ct, c);
		Multiply(out33, ct, c);
		out33.Display();
		tcout << tendl;


		//  c * c^t =
		//   [  107   137 ]
		//   [ -137   179 ]

		tcout << _T( "MultiplyTranspose(out, c, c);\n" );
		Zero(out22);
		MultiplyTranspose(out22, c, c);
		out22.Display();
		tcout << tendl;

		tcout << _T( " MultiplyTranspose(out,  c);\n" );
		Zero(out22);
		MultiplyTranspose(out22, c);
		out22.Display();
		tcout << tendl;

		tcout << _T( "Transpose(ct, c) and Multiply(out, c, ct);\n" );
		Zero(out22);
		Multiply(out22, c, ct);
		out22.Display();
		tcout << tendl;


		buf.Display();
	}
	

	return 0;
}