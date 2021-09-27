#include	<chrono>
#include	<crtdbg.h>


#include	<oreore/mathlib/GraphicsMath.h>

#include	<matrixlib/MatrixLib.h>




int main()
{
	_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	{
		DynamicMatrix<double> rgbimage( 16384, 16384 );

		std::chrono::system_clock::time_point  start, end; // 型は auto で可
		start = std::chrono::system_clock::now(); // 計測開始時間

		for( int y=0; y<rgbimage.numRows(); ++y )
		{
			for( int x=0; x<rgbimage.numCols(); ++x )
			{
				double& pixel = rgbimage( y, x );
				pixel = 255.0;
				//Vec4f& pixel = rgbimage( y, x);
				//InitVec( pixel, 255.0f, 0.5f, 0.0f, 1.0f );
			}
		}

		end = std::chrono::system_clock::now();  // 計測終了時間
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( end-start ).count(); //処理に要した時間をミリ秒に変換
		tcout << "time elapsed: " << elapsed << "[ms].\n";
	
	}


	return 0;
}