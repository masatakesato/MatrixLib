#include	<oreore/mathlib/Random.h>

#include	"MovingLeastSquares.h"


#define NUM_X 3
#define	NUM_Y 3





int main( int argc, char** argv )
{

	Controller controller(9);

	StaticMeshGrid<10, 10> mesh;


	Matrix<float, NUM_X*NUM_X, 2> controlPoints;

	for( int y=0; y<NUM_Y; ++y )
	{
		for( int x=0; x<NUM_X; ++x )
		{
			//controller.MoveTo( y*NUM_X+x, float(x), float(y) );
			controlPoints( y*NUM_X+x, 0 ) = float(x);
			controlPoints( y*NUM_X+x, 1 ) = float(y);
		}
	}



	Matrix<float, NUM_X*NUM_X, 2> targetPoints;

	for( int y=0; y<NUM_Y; ++y )
		for( int x=0; x<NUM_X; ++x )
		{
			targetPoints( y*NUM_X+x, 0 ) = float(x);
			targetPoints( y*NUM_X+x, 1 ) = float(y);
		}

	targetPoints(NUM_X + 1, 0) = 1.5f;
	targetPoints(NUM_X + 1, 1) = 1.0f;

//	InitVec( targetPoints[ NUM_X + 1 ], 1.5f, 1.0f );


	MovingLeastSquares mls;

	//Matrix<float> w(1, controlPoints.numRows() );
	//fArray A(9);
	//OreOreLib::Array<Vec2f> p;
	OreOreLib::Array<Vec2f> vs(1);
	Vec2f v2;
	InitVec( vs[0], 1.0f, 1.0f );
	float alpha = 0.01f;

	//for( int i=0; i<1200*1200; ++i )
	{
		mls.PrecomputeAffineAs( controlPoints, vs, alpha );

		//mls.PrecomputeAffineA( w, A, controlPoints, v, a );
		//mls.TransformAffine( targetPoints, w, A ).Display();
	}


//	Vec2f v, v2;
//	InitVec(v, 1.0f, 1.0f);
//	mls.Transform( v2, v );

	//tcout << tendl;

	//for( double r=-1.0; r<=1.0; r+=0.1 )
	//{
	//	tcout << r << ": " << r*r * log(abs(r)) << tendl;
	//}

	//tcout << log(r) << tendl;
	return 0;
}