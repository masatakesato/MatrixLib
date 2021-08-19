#ifndef MOVING_LEAST_SQUARES_H
#define	MOVING_LEAST_SQUARES_H


#include	<oreore/container/StaticArray.h>
//#include	<oreore/mathlib/GraphicsMath.h>//#include	<oreore/mathlib/MathLib.h>


#include	<matrixlib/MatrixLib.h>
#include	<matrixlib/VectorOperations.h>
#include	<matrixlib/LU.h>


#include	"IDeformer.h"
#include	"IMeshGrid.h"




class Controller

{
public:

	Controller()
	{
	}


	Controller( int numhandles )
	{
		m_Origin.Init( numhandles );
	}


	~Controller()
	{
		m_Origin.Release();
	}


	int NumHandles() const
	{
		return m_Origin.Length();
	}


	bool AddHandle()
	{
		m_Origin.AddToTail();
	}


	bool RemoveHandle( int i=0 )
	{
		m_Origin.Remove(i);
	}


	Vec2f HandleOrigin( int i ) const
	{
		return m_Origin[i];
	}


	void SetHandleOrigin( int i, float x, float y )
	{
		InitVec( m_Origin[i], x, y );
	}


	void OffsetHandleOrigin( int i, float dx, float dy )
	{
		m_Origin[i](0) += dx;
		m_Origin[i](1) += dy;
	}


	Vec2f HandleTarget( int i ) const
	{
		return m_Target[i];
	}

	void SetHandleTarget( int i, float x, float y )
	{
		InitVec( m_Origin[i], x, y );
	}


	void OffsetHandleTarget( int i, float dx, float dy )
	{
		m_Target[i](0) += dx;
		m_Target[i](1) += dy;
	}


	const OreOreLib::Array<Vec2f> HandleOrigins() const
	{
		return m_Origin;
	}


	const OreOreLib::Array<Vec2f> HandleTargets() const
	{
		return m_Origin;
	}



private:

	OreOreLib::Array<Vec2f>	m_Origin;
	OreOreLib::Array<Vec2f>	m_Target;


//	拡張性を考えると配列(というかリスト)の方がよいが、
//	計算の利便性考えると行列の方がよい。
//		↓
//	計算が必要な際に行列にデータコピーする
};





class MeshGrid : public IMeshGrid
{
public:

	MeshGrid()
		: m_DimX(0)
		, m_DimY(0)
	{
	}


	MeshGrid( int xdim, int ydim )
		: m_DimX(xdim)
		, m_DimY(ydim)
	{
		m_Vertices.Init( m_DimX * m_DimY, 2 );
	}


	~MeshGrid()
	{
		m_Vertices.Release();
	}


	const Vec2f& At( int x, int y )
	{
		return (Vec2f&)m_Vertices( y*m_DimY + x );
	}


	int DimX() const
	{
		return m_DimX;
	}


	int DimY() const
	{
		return m_DimY;
	}


	int numVertices() const
	{
		return m_DimX * m_DimY;
	}


	const DynamicMatrix<float>& V() const
	{
		return m_Vertices;
	}



private:

	int m_DimX, m_DimY;
	DynamicMatrix<float>	m_Vertices;	// num of vertices * dim(=2)


};



template< size_t DimX, size_t DimY >
class StaticMeshGrid : public IMeshGrid
{
public:

	StaticMeshGrid()
	{
	}


	~StaticMeshGrid()
	{
	}


	const Vec2f& At( int x, int y )
	{
		return m_Vertices[ y*DimY + x ];
	}


private:

	OreOreLib::StaticArray<Vec2f, DimX * DimY> m_Vertices;


};




class MovingLeastSquares
{

public:

	MovingLeastSquares();
	virtual ~MovingLeastSquares();

	void Init( int numControlPoints, int numVertices );



	void PrecomputeAffineA( IMatrix<float>& w, IMatrix<float>& A, const IMatrix<float>& p, const Vec2f& v, float alpha );
	Vec2f TransformAffine( const IMatrix<float>& q, IMatrix<float>& w, const OreOreLib::Array<float>& A );


	void PrecomputeAffineAs( const IMatrix<float>& p, const OreOreLib::Array<Vec2f>& v, float alpha );


private:

	void PrecomputeWeights( IMatrix<float>& w, const IMatrix<float>& p, const Vec2f& v, float alpha );
	void ComputeWCentroids( Vec2f& p_star, const IMatrix<float>& p, const IMatrix<float>& w );
	void ComputeHats( IMatrix<float>& p_hat, const IMatrix<float>& p, const Vec2f& p_star );



	int m_numControlPoints;
	int m_numVertices;



};


//Aとwは、補完したい点vごとに必要になる

//メッシュグリッドの頂点数は頻繁に変化するか？ -> No
//制御点の数は頻繁に変化するか？ -> Yes


//Aとwは、メッシュグリッド頂点毎に必要.



class MLSDeformer : public IDeformer
{
public:

	MLSDeformer( MeshGrid& mesh, Controller& ctrl )
		: m_refMesh( mesh )
		, m_refController( ctrl )
	{
	}


	~MLSDeformer(){}


	void Update();



private:

	MeshGrid&	m_refMesh;
	Controller& m_refController;

	// Intermediate buffers for Affine
	DynamicMatrix<float>	m_W,		// num of p * num of v
							m_A,		// num of p * num of v
							m_Pstar,	// num of v * dim(=2)
							m_Ph_x,		// num of p * num of v
							m_Ph_y;		// num of p * num of v 

	MovingLeastSquares		m_MLSSolver;

};





#endif // !MOVING_LEAST_SQUARES_H

