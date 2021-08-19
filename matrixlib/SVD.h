#ifndef SVD_H
#define	SVD_H


#include	"Matrix.h"



// 固有値分解(ヤコビ法)
// a: 対称行列. 対角要素に"固有値"を格納する行列
// v: 固有ベクトル
// eps: 収束誤差
// iter_max: 反復回数の最大値
// i行j列

// out : V, Lamda, in : A で引数整理
/*
template< typename T >
int Eigen_JacobiMethod( Matrix<T>& a, Matrix<T>& v, double eps=1.0e-8, int iter_max=100)
{
	if( !a.IsSquare() )	return -1;	//正方行列だけ扱う

	const int	N = a.numCols();

    double *bim, *bjm;
    double bii, bij, bjj, bji;
 
	bim = new double[N];
    bjm = new double[N];
 
	Identity( v );
 
    int cnt = 0;
    for(;;)
	{
        int i, j;
 
        T x = 0.0;
        for(int ia = 0; ia < N; ++ia)
		{
            for(int ja = 0; ja < N; ++ja)
			{
                int idx = ia*N+ja;
				if( ia != ja && fabs( *a.ptr(ia, ja) ) > x )
				{
                    i = ia;
                    j = ja;
					x = fabs( *a.ptr(ia, ja) );
                }
            }
        }
 
        double aii = *a.ptr(i, i);
        double ajj = *a.ptr(j, j);
        double aij = *a.ptr(i, j);
		
        double alpha, beta;
        alpha = (aii-ajj)/2.0;
        beta  = sqrt(alpha*alpha+aij*aij);
 
        double st, ct;
        ct = sqrt((1.0+fabs(alpha)/beta)/2.0);    // sinθ
        st = (((aii-ajj) >= 0.0) ? 1.0 : -1.0)*aij/(2.0*beta*ct);    // cosθ
 
        // A = PAPの計算
        for(int m = 0; m < N; ++m)
		{
            if(m == i || m == j) continue;
 
            double aim = *a.ptr(i, m);
            double ajm = *a.ptr(j, m);
 
            bim[m] =  aim*ct + ajm*st;
            bjm[m] = -aim*st + ajm*ct;
        }
 
        bii = aii*ct*ct + 2.0*aij*ct*st + ajj*st*st;
        bij = 0.0;
 
        bjj = aii*st*st - 2.0*aij*ct*st + ajj*ct*ct;
        bji = 0.0;
 
        for(int m = 0; m < N; ++m)
		{
            *a.ptr(i, m) = *a.ptr(m, i) = bim[m];
            *a.ptr(j, m) = *a.ptr(m, j) = bjm[m];
        }

        *a.ptr(i, i) = bii;
        *a.ptr(i, j) = bij;
        *a.ptr(j, j) = bjj;
        *a.ptr(j, i) = bji;
 
        // V = PVの計算
        for(int m = 0; m < N; ++m)
		{
            double vmi = *v.ptr(m, i);
            double vmj = *v.ptr(m, j);
 
            bim[m] =  vmi*ct + vmj*st;
            bjm[m] = -vmi*st + vmj*ct;
        }
        for(int m = 0; m < N; ++m)
		{
            *v.ptr(m, i) = bim[m];
            *v.ptr(m, j) = bjm[m];
        }
 
        double e = 0.0;
        for(int ja = 0; ja < N; ++ja)
		{
            for(int ia = 0; ia < N; ++ia)
			{
                if(ia != ja)
				{
                    e += fabs( *a.ptr(ja, ia) );
                }
            }
        }
        if(e < eps) break;
 
        cnt++;
        if(cnt > iter_max) break;
    }// end of while
 
    delete [] bim;
    delete [] bjm;
 
    return cnt;
}
*/



// Singular Value Decomposition
// U
template< typename T >
void SVD( Matrix<T> &U, Matrix<T> &S, Matrix<T> &V, const Matrix<T> &A, double eps=1.0e-8 )
{

	Matrix<T>	At, B;

	// B = A * A_transpose計算
	Transpose( At, A );
	Multiply( B, A, At );

	// Bの固有ベクトル Uを計算する
	Eigen_JacobiMethod( B, U, eps );


}


#endif	// SVD_H //