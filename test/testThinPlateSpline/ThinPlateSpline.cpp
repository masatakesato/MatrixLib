#include	"ThinPlateSpline.h"


// https://khanhha.github.io/posts/Thin-Plate-Splines-Warping/

// https://vision.cornell.edu/se3/wp-content/uploads/2014/09/cvpr05b.pdf


static float tps_kernel( const Vec2f& p, const Vec2f& q )
{
	float r2 = Max( DistanceSqrd( p, q ), FLT_EPSILON );
	return r2 * log( sqrt(r2) );
}


template< typename T >
static float tps_kernel( const T p[], const T q[], const int dim )
{
	T r2 = 0;
	for( int i=0; i<dim; ++i )
		r2 += pow( p[i] - q[i], 2 );
	
	return T( r2 * log( sqrt(Max(r2, (T)1.0e-9))) );
}




ThinPlateSpline::ThinPlateSpline()
{
	
}
	

ThinPlateSpline::~ThinPlateSpline()
{

}



void ThinPlateSpline::Init( int numDim, int numControlPoints )
{
	assert(numDim>0);
	assert(numControlPoints>0);

	m_Dim				= numDim;
	m_NumControlPoints	= numControlPoints;
	m_MatSize			= numControlPoints + 1 + numDim;

	m_A.Init( m_MatSize, m_MatSize );
	m_x.Init( m_MatSize, m_Dim );
	m_b.Init( m_MatSize, m_Dim );

	m_K.Init( m_A, 0, 0, m_NumControlPoints, m_NumControlPoints );	// Reference to K part of TPS matrix
	m_P.Init( m_A, 0, m_NumControlPoints+1, m_NumControlPoints, m_Dim );// Reference to P part of TPS matrix (except first column filled with 1)
	m_Pt.Init( m_A, m_NumControlPoints+1, 0, m_Dim, m_NumControlPoints );// Reference to P^t part of TPS matrix (except first column filled with 1)

	m_w.Init( m_x, 0, 0, m_NumControlPoints, m_Dim );// Reference to w part of x
	m_a.Init( m_x, m_NumControlPoints, 0, m_Dim+1, m_Dim );// Reference to a part of x

	m_U.Init( 1, m_NumControlPoints );


	// fill-in fixed values to A
	for( int i=0; i<m_NumControlPoints; ++i )
	{
		m_A(m_NumControlPoints, i) = 1;
		m_A(i, m_NumControlPoints) = 1;
	}

	m_Solver.Init( m_A, numDim );
}



void ThinPlateSpline::Release()
{
	m_K.Release();
	m_P.Release();
	m_Pt.Release();

	m_A.Release();
	m_x.Release();
	m_b.Release();

	m_U.Release();

	//m_Solver.Release();
}



void ThinPlateSpline::Update( OreOreLib::Array<Vec2f>& controlPoints )
{	
	// K
	for( int i=0; i<m_K.numRows(); ++i )
	{
		// non-diagonal elements
		for( int j=0; j<i; ++j )
		{
			auto u = tps_kernel( controlPoints[i], controlPoints[j] );
			m_K(i, j) = u;
			m_K(j, i) = u;
		}

		// diagonal elements
		m_K(i, i) = 0.0f;//or regularizationParameter;//
	}
	
	// P and Pt
	for( int i=0; i<m_NumControlPoints; ++i )
	{
		// P
		m_P(i, 0) = controlPoints[i].x;
		m_P(i, 1) = controlPoints[i].y;

		// P^t
		m_Pt(0, i) = controlPoints[i].x;
		m_Pt(1, i) = controlPoints[i].y;
	}

	m_A.Display();

	m_Solver.SetCoeff( m_A );
}



void ThinPlateSpline::Solve( OreOreLib::Array<Vec2f>& targetPoints )
{
	for( int i=0; i<Min( targetPoints.Length(), m_NumControlPoints ); ++i )
	{
		m_b(i, 0) = targetPoints[i].x;
		m_b(i, 1) = targetPoints[i].y;
	}

//	m_b.Display();

	m_Solver.Solve( m_x, m_b );

//	m_x.Display();
//	m_x.Display();
//	m_w.Display();
//	m_a.Display();
}



void ThinPlateSpline::Transform( Vec2f& out, const Vec2f& in )
{
	// x only
	float (&vec)[2] = out.xy;
	const float (&invec)[2] = in.xy;


	for( int i=0; i<m_U.numElms(); ++i )
		m_U(i) = tps_kernel( m_P.ptr(i, 0), invec, m_Dim );


	for( int d=0; d<m_Dim; ++d )
	{
		// Linear part
		vec[d] = m_a(0, d);
		for( int j=0; j<m_a.numRows()-1; ++j )
			vec[d] += invec[j] * m_a(j+1, d);

		// Kernel
		for( int j=0; j<m_NumControlPoints; ++j )
			vec[d] += m_w(j, d) * m_U(j);
	}

//	tcout << tendl << out.x << ", " << out.y << tendl;

}