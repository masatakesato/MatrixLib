#include	"MovingLeastSquares.h"



MovingLeastSquares::MovingLeastSquares()
	: m_numControlPoints(0)
	, m_numVertices(0)
{
	
}
	

MovingLeastSquares::~MovingLeastSquares()
{

}




void MovingLeastSquares::Init( int numControlPoints, int numVertices )
{
	m_numControlPoints = numControlPoints;
	m_numVertices	= numVertices;
}



/*

Moving Least Squares. Affine
v: 任意の座標
p: 初期の制御点座標
q: 移動後の制御点座標

wi: v～pi間の重み。  1.0 / | p[i] - v | ^(2*α)
p* = Σ wi*pi / Σ wi
q* = Σ wi*qi / Σ wi

pi^ = pi - p*
qi^ = qi - q*

Aj = (v-p*) * (  ??? )^-1 * ???2


[ p0x, p0y ]
[ p1x, p1y ]
[ p2x, p2y ]
[ p3x, p3y ]

//============= ??? ==============//
??? = 
[ p0x ] * w0 * [ p0x, p0y ] = 
[ p0y ]

Σ  wi * [  pix*pix, pix*piy ]
           [  piy*pix, piy*piy ]


2x2_inv = 
[ a22, -a12 ]
[ -a21, a11 ]
/ (a11*a22 -a12*a21)

//============== ???2 =============//


mat1x2 * mat2x2 * mat2x1


*/



void MovingLeastSquares::PrecomputeWeights( IMatrix<float>& w, const IMatrix<float>& p, const Vec2f& v, float alpha )
{

	for( int i=0; i<p.numRows(); ++i )
	{
		auto norm = DotRR( p, v, i, 0 );
		w(i) = 1 / pow( Max( norm, (std::numeric_limits<float>::min)() ), alpha );

		//w(i, j) = 1 / pow( Distance( p[i], v(j) ), 2*alpha );
	}
}



void MovingLeastSquares::ComputeWCentroids( Vec2f& p_star, const IMatrix<float>& p, const IMatrix<float>& w )
{
	TransposeMultiply( p_star, w, p );	// Σ wi * pi
	Scale( p_star, 1 / SumCol( w ) );	// Σ wi * pi / Σ wi
}



void MovingLeastSquares::ComputeHats( IMatrix<float>& p_hats, const IMatrix<float>& p, const Vec2f& p_star )
{
	for( int i=0; i<p.numRows(); ++i )
	{
		for( int j=0; j<p.numCols(); ++j )
		{
			p_hats(i, j) = p(i, j) - p_star(j);
		}
	}
}



void MovingLeastSquares::PrecomputeAffineA( IMatrix<float>& w, IMatrix<float>& A, const IMatrix<float>& p, const Vec2f& v, float alpha )
{
	PrecomputeWeights( w, p, v, alpha );

	Vec2f p_star;
	ComputeWCentroids( p_star, p, w );// Can be precomputed from p and v

	DynamicMatrix<float> p_hats( p.numRows(), p.numCols() );// = p;
	ComputeHats( p_hats, p, p_star );// Can be precomputed from p and v

	Matrix<float, 2, 2> ptwp, ptwp_inv;
	float a11=0, a12=0, a22=0;

	//Σ  wi * [ pix*pix, pix*piy ]   -> a11= wi * pix*pix; a12 = wi * pix*piy; a22 = wi * piy*piy;
	//        [ piy*pix, piy*piy ]
	for( int i=0; i<p_hats.numRows(); ++i )
	{
		float wi = w(i);
		a11 += wi * pow( p_hats(i, 0), 2 );
		a12 += wi * p_hats(i, 0) * p_hats(i, 1);
		a22 += wi * pow( p_hats(i, 1), 2 );
	}

	// ptwp_inv
	//2x2_inv = 
	// [ a22, -a12 ]  / (a11*a22 -a12*a21)
	// [ -a21, a11 ]
	auto det_inv = 1 / ( a11*a22 - a12*a12 );

	ptwp_inv(0, 0) = a22 * det_inv;		ptwp_inv(0, 1) = -a12 * det_inv;
	ptwp_inv(1, 0) = -a12 * det_inv;	ptwp_inv(1, 1) = a11 * det_inv;

	Vec2f v_hat;
	Subtract( v_hat, v, p_star );
 
	// (v - p*) * (Σ p_hat^t * w * p^hat )^-1
	Vec2f tmp;
	tmp(0) = v_hat(0) * ptwp_inv(0, 0) + v_hat(1) * ptwp_inv(1, 0);
	tmp(1) = v_hat(0) * ptwp_inv(0, 1) + v_hat(1) * ptwp_inv(1, 1);

	// * ( wj * pj^t )
	for( int j=0; j<A.numElms(); ++j )
	{
		//wp.x = w(j) * p_hat[j].x;
		//wp.y = w(j) * p_hat[j].y;

		A(j) = tmp(0) * w(j) * p_hats(j, 0)
				+ tmp(1) * w(j) * p_hats(j, 1);
	}

}





void MovingLeastSquares::PrecomputeAffineAs( const IMatrix<float>& p, const OreOreLib::Array<Vec2f>& v, float alpha )
{
	// Precompute w
	DynamicMatrix<float> As( p.numRows(), v.Length() );
	DynamicMatrix<float> Ws( p.numRows(), v.Length() );

	//for( int i=0; i<v.Length(); ++i )
	//{
	//	MatrixView<float> A( As, 0, i, p.numRows(), 1 );
	//	MatrixView<float> w( Ws, 0, i, p.numRows(), 1 );
	//	
	//	PrecomputeAffineA( w, A, p, v[i], alpha );
	//}
}





Vec2f MovingLeastSquares::TransformAffine( const IMatrix<float>& q, IMatrix<float>& w, const OreOreLib::Array<float>& A )
{
	Vec2f q_star;
	ComputeWCentroids( q_star, q, w );
	
	DynamicMatrix<float> q_hats( q.numRows(), q.numCols() );
	ComputeHats( q_hats, q, q_star );

	Vec2f Fa = q_star;

	// Accumulate Σ (A(j) * q_hat(j))
	for( int j=0; j<q.numRows(); ++j )
	{

		Fa(0) += q_hats(j, 0) * A[j];
		Fa(1) += q_hats(j, 1) * A[j];
		//A[j] * 
		//AddScaled( Fa, q_hats(j], A[j] );
	}

	return Fa;
}




void MLSDeformer::Update()
{
	MatrixView<float> w, A, p;

//	p.Init( m_refController.Origin().begin(), m_refController.NumHandles(), 2 );

	for( int i=0; i<m_refMesh.numVertices(); ++i )
	{
		w.Init( m_W, 0, i, m_refController.NumHandles(), 1 );
		A.Init( m_A, 0, i, m_refController.NumHandles(), 1 );

//		m_MLSSolver.PrecomputeAffineA( w, A, p, v, alpha )

	}

}