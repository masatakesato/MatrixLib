#ifndef VARIANCE_COVARIANCE_CORRELATION_H
#define VARIANCE_COVARIANCE_CORRELATION_H


#include "MatrixLib.h"


//extern double Variance_Covariance(double *a, double *b, unsigned int num);


// https://www.ouj.ac.jp/mijika/tokei/xml/k3_06022.xml


// Calculate Covariance Matrix.
// in: sample vector data with row-major order
// 
template< typename T >
void Covariance( IMatrix<T>& out, const IMatrix<T>& in )
{
	T num_inv = 1.0 / in.numRows();

	for( int i=0; i<out.numRows(); i++ )
	{
		// Calculate average of i-th dimensional emelents from in.
		T avg_i = 0;
		for( int k=0; k<in.numRows(); k++ )	avg_i += in( k, i );
		avg_i *= num_inv;

		for( int j=i; j<out.numCols(); j++ )
		{
			// Calculate average of j-th dimensional emelents from in.
			T avg_j = 0;
			if( i==j )
			{ 
				avg_j = avg_i;
			}
			else
			{
				for( int k=0; k<in.numRows(); k++ )	avg_j += in( k, j );
				avg_j *= num_inv;
			}

			// Calculate covariance between dimension i and j.
			T& cov = out( i, j );
			cov = 0;
			for( int k=0; k<in.numRows(); k++ )	cov += ( in( k, j ) - avg_j ) * ( in( k, i ) - avg_i );
			cov *= num_inv;

			// Copy covariance to diagonal element
			out( j, i ) = cov;

		}// end of j loop
	}// end of i loop

}



/*
// 相関
double Correlation(double *a, double *b, unsigned int num)
{

return 0;//Variance_Covariance(a,b,num) / sqrt( Variance_Covariance(a,a,num) * Variance_Covariance(b,b,num) );

}
*/





/*
// Covariance of vector a and b
double Variance_Covariance(double *a, double *b, unsigned int num){// aとbの共分散を求める。aとbが同じ配列の時は分散

unsigned int i;
double num_inv = 1.0/num,
v = 0,		//variance/covariance
a_ave = 0,	//aの平均
b_ave = 0;	//bの平均

//************** aとbの平均を計算 ***************
for(i=0; i<num; i++){
a_ave += *(a+i);  b_ave += *(b+i);
}
a_ave *= num_inv;  b_ave *= num_inv;

//************** (ai - a_ave)(bi - b_ave) **************
for(i=0; i<num; i++) v += ( *(a+i) - a_ave ) * ( *(b+i) - b_ave );

printf("%lf\n", v * num_inv);

return v *= num_inv;
}



*/



#endif /* VARIANCE_COVARIANCE_CORRELATION_H */