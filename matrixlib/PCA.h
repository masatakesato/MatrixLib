#ifndef PCA_H
#define	PCA_H


#include	"MatrixLib.h"



// 主成分分析クラス
template< typename T>
class PCA
{
private:

	DynamicMatrix<T>	C;
	DynamicMatrix<T>	A1;
	DynamicMatrix<T>	A2;
	DynamicMatrix<T>	X1;
	DynamicMatrix<T>	X2;


	void Allocate();
	void Release();

	int Jacobi(	int n, int ct, double eps, IMatrix<T>& A, IMatrix<T>& A1, IMatrix<T>& A2, IMatrix<T>& X1, IMatrix<T>& X2 );


public:

	PCA(){};
	~PCA(){};

	int Principal( int p, int n, IMatrix<T>& x, double r[], IMatrix<T>& a, double eps, int ct );


	void Execute();
};






#endif // !PCA_H
