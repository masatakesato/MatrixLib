#ifndef THIN_PLATE_SPLINE_H
#define	THIN_PLATE_SPLINE_H


#include	<oreore/mathlib/GraphicsMath.h>

#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/LU.h>



class ThinPlateSpline
{

public:

	ThinPlateSpline();
	virtual ~ThinPlateSpline();

	void Init( int numDim, int numControlPoints );
	void Release();// TODO: Implement

	void Update( OreOreLib::Array<Vec2f>& controlPoints );

	void Solve( OreOreLib::Array<Vec2f>& targetPoints );

	void Transform( Vec2f& out, const Vec2f& in );


protected:

	DynamicMatrix<float>	m_A, m_x, m_b;	// L, P, Pt

	MatrixView<float>	m_K, m_P, m_Pt,
						m_w, m_a;

	PivotLU_Solver<float>	m_Solver;

	int	m_Dim;
	int	m_NumControlPoints;
	int	m_MatSize;

	DynamicMatrix<float>	m_U;



private:


};



#endif // !THIN_PLATE_SPLINE_H

