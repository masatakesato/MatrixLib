//#include "VarianceCovariance_Correlation.h"
//
//
//#include <stdio.h>
//#include <math.h>
//
//
//
//double Variance_Covariance(double *a, double *b, unsigned int num){// aとbの共分散を求める。aとbが同じ配列の時は分散
//
//	unsigned int i;
//	double num_inv = 1.0/num,
//		   v = 0,		//variance/covariance
//		   a_ave = 0,	//aの平均
//		   b_ave = 0;	//bの平均
//
//	//************** aとbの平均を計算 ***************
//	for(i=0; i<num; i++){
//		a_ave += *(a+i);  b_ave += *(b+i);
//	}
//	a_ave *= num_inv;  b_ave *= num_inv;
//
//	//************** (ai - a_ave)(bi - b_ave) **************
//	for(i=0; i<num; i++) v += ( *(a+i) - a_ave ) * ( *(b+i) - b_ave );
//
//	printf("%lf\n", v * num_inv);
//
// return v *= num_inv;
//}
//
//
//void Variance_Covariance_Matrix(Matrix *dest, Matrix *a){
//
//	double ave_i=0, ave_j=0, num_inv = 1.0/a->Column;
//	unsigned int i,j,k;
//
//	for(i=0; i<dest->Row; i++){// i行
//
//		/******************* i行目の平均 *******************/
//		ave_i = 0;
//		for(k=0; k<a->Column; k++) ave_i += *(a->Elements + i*a->Column +k);//
//		ave_i *= num_inv;
//
//		for(j=i; j<dest->Column; j++){// i列
//
//			/*******************j行目の平均 *******************/
//			ave_j = 0;
//			for(k=0; k<a->Column; k++) ave_j += *(a->Elements + j*a->Column +k);//j行目の平均
//			ave_j *= num_inv;
//
//			/*(i行目-i_ave)*(j行目-j_ave)*/
//			for(k=0; k<a->Column; k++) 
//				*(dest->Elements + i*dest->Column +j) += ( *(a->Elements + i*a->Column +k) - ave_i ) * ( *(a->Elements + j*a->Column +k) - ave_j);
//			
//			*(dest->Elements + i*dest->Column +j) *=num_inv;
//			*(dest->Elements + j*dest->Column +i) = *(dest->Elements + i*dest->Column +j);
//
//		
//		}// end of j loop
//	}// end of i loop
//
//}
//
//
//
//
///*
//double Correlation(double *a, double *b, unsigned int num){
//
//	return 0;//Variance_Covariance(a,b,num) / sqrt( Variance_Covariance(a,a,num) * Variance_Covariance(b,b,num) );
//
//}
//*/