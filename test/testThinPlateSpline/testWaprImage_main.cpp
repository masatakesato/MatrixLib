#include	<FreeImage.h>

#include	"ThinPlateSpline.h"



#define NUM_X 3
#define	NUM_Y 3




int main()
{

	FreeImage_Initialise();

	//====================== Read Image ======================//

	FREE_IMAGE_FORMAT fif	= FreeImage_GetFileType("V:/sampleimages/image.png", 0);
	FIBITMAP* bitmap = FreeImage_Load( fif, "V:/sampleimages/image.png" );

	unsigned width	= FreeImage_GetWidth( bitmap );
	unsigned height	= FreeImage_GetHeight( bitmap );

	

	//================ Initialize Output Image ===============//

	unsigned bpp = FreeImage_GetBPP( bitmap );
	FIBITMAP* outbitmap  = FreeImage_Allocate( width, height, bpp );



	//====================== Warp Image ======================//


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

	InitVec( targetPoints[ NUM_X + 1 ], 1.0f, 1.5f );


	ThinPlateSpline tps;

	tps.Init(2, controlPoints.Length() );
	tps.Update( controlPoints );
	tps.Solve( targetPoints );


	RGBQUAD pixcolor;

	for( unsigned y=0; y<height; ++y )
	{
		for( unsigned x=0; x<width; ++x )
		{

			Vec2f v, v2, diff;
			InitVec(v, 2.0f * float(x)/float(width), 2.0f * float(y)/float(width) );
			tps.Transform( v2, v );

			Subtract( diff, v, v2 );
			Add( v, diff );
			auto x_ = unsigned(v2.x*width);
			auto y_	= unsigned(v2.y*height);
			FreeImage_GetPixelColor( bitmap, Clamp(unsigned(v.x*width/2), unsigned(0), width-1), Clamp(unsigned(v.y*height/2), unsigned(0), height-1), &pixcolor );
			FreeImage_SetPixelColor( outbitmap, x, y, &pixcolor );
		}
	}


	//===================== Output ===========================//
	FreeImage_Save( fif, outbitmap, "./warpedimage.png" );



	//==================== Release =======================//
	FreeImage_Unload( bitmap );
	FreeImage_Unload( outbitmap );

	FreeImage_DeInitialise(); 

	return 0;
}