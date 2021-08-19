#include	<oreore/mathlib/Random.h>

#include	"ThinPlateSpline.h"


#define NUM_X 3
#define	NUM_Y 3



int main( int argc, char** argv )
{
	OreOreLib::Array<Vec2f> controlPoints(NUM_X * NUM_Y);

	for( int y=0; y<NUM_Y; ++y )
		for( int x=0; x<NUM_X; ++x )
			InitVec( controlPoints[y*NUM_X+x], float(x), float(y) );


	OreOreLib::Array<Vec2f> targetPoints(NUM_X * NUM_Y);

	for( int y=0; y<NUM_Y; ++y )
		for( int x=0; x<NUM_X; ++x )
		{
			const Vec2f& o = controlPoints[y*NUM_X+x];
			InitVec( targetPoints[y*NUM_X+x], float(x), float(y) );
		}

	InitVec( targetPoints[ NUM_X + 1 ], 1.5f, 0.5f );


	ThinPlateSpline tps;

	tps.Init(2, controlPoints.Length() );


	tps.Update( controlPoints );

	tps.Solve( targetPoints );


	Vec2f v, v2;
	InitVec(v, 1.0f, 1.0f);
	tps.Transform( v2, v );

	//tcout << tendl;

	//for( double r=-1.0; r<=1.0; r+=0.1 )
	//{
	//	tcout << r << ": " << r*r * log(abs(r)) << tendl;
	//}

	//tcout << log(r) << tendl;
	return 0;
}