#ifndef SOLVER_H
#define	SOLVER_H



#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP


#include	"MatrixLib.h"


template< typename T >
class Solver
{
public:

	Solver(){}
	virtual ~Solver(){}


	virtual void Init( const IMatrix<T>& A, int batchSize=1 ){}
	virtual bool SetCoeff( const IMatrix<T>& A )=0;// Updates coefficnet matrix without reallocating
	virtual bool SetBatchSize( int batchSize )=0;

	virtual void Solve( IMatrix<T>& x, const IMatrix<T>& b )=0;
};




#endif // !SOLVER_H
