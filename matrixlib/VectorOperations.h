#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include	<oreore/common/Utility.h>
#include	<oreore/mathlib/MathLib.h>

#include	"StaticMatrix.h"





template< typename T >
using Vec2		= StaticMatrix<T, 1, 2>;
//template< typename T > union	Vec2;

using Vec2uc	= Vec2<unsigned char>;
using Vec2s		= Vec2<short>;
using Vec2us	= Vec2<unsigned short>;
using Vec2i		= Vec2<int>;
using Vec2ui	= Vec2<unsigned int>;
using Vec2f		= Vec2<float>;
using Vec2d		= Vec2<double>;


template< typename T >
using Vec3		= StaticMatrix<T, 1, 3>;
//template< typename T > union	Vec3;

typedef Vec3<unsigned char>		Vec3uc;
typedef Vec3<short>				Vec3s;
typedef Vec3<unsigned short>	Vec3us;
typedef Vec3<int>				Vec3i;
typedef Vec3<unsigned int>		Vec3ui;
typedef Vec3<float>				Vec3f;
typedef Vec3<double>			Vec3d;


template< typename T >
using Vec4 = StaticMatrix<T, 1, 4>;
//template< typename T > union	Vec4;

typedef Vec4<unsigned char>		Vec4uc;
typedef Vec4<short>				Vec4s;
typedef Vec4<unsigned short>	Vec4us;
typedef Vec4<int>				Vec4i;
typedef Vec4<unsigned int>		Vec4ui;
typedef Vec4<float>				Vec4f;
typedef Vec4<double>			Vec4d;


template< typename T > union	Mat4;

typedef Mat4<float>				Mat4f;
typedef Mat4<double>			Mat4d;


template< typename T > union	Quaternion;

typedef	Quaternion<float>		Quatf;
typedef	Quaternion<double>		Quatd;




//##############################################################################//
//										Scalar									//
//##############################################################################//


union ieee754
{
	struct
	{
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign     : 1;
	};

	float            _f32;
	unsigned int     _u32;
};



//##############################################################################//
//									2D Vector									//
//##############################################################################//

// 2次元ベクトル共用体
//template< typename T >
//union Vec2
//{
//
//	struct { T x, y; };
//	struct { T u, v; };
//	T	xy[2];
//	T	uv[2];
//
//	Vec2(): x( 0 ), y( 0 ) {}
//	Vec2( T x_, T y_ ): x( x_ ), y( y_ ) {}
//
//	//=========== experimental implementation. 2018.10.14 ============//
//	// Copy constructor
//	Vec2( const Vec2& obj )
//	{
//		x = obj(0);
//		y = obj(1);
//	}
//
//	// Copy constructor
//	~Vec2()
//	{
//	}
//
//	// Move constructor
//	Vec2( Vec2&& obj )
//	{
//		x = obj(0);
//		y = obj(1);
//	}
//
//	// Copy assignment operator
//	Vec2& operator=( const Vec2& obj )
//	{
//		if( this != &obj )
//		{
//			x = obj(0);
//			y = obj(1);
//		}
//		return *this;
//	}
//
//	// Move assignment opertor =
//	Vec2& operator=( Vec2&& obj )
//	{
//		if( this != &obj )
//		{
//			x = obj(0);
//			y = obj(1);
//		}
//		return *this;
//	}
//
//	T& operator()( int i ) &
//	{
//		return xy[i];
//	}
//
//	// Subscription operator. ( called by following cases: "T& a = Matrix<T>(10,10)(n)", "auto&& a = Matrix<T>(20,4)(n)" )
//	T operator()( int i ) const&&
//	{
//		return std::move(xy[i]);
//	}
//
//};



// Init Vector
template< typename T >
void InitVec( Vec2<T>& inout, T x, T y )
{
	inout(0)	= x;
	inout(1)	= y;
}


// Init Vector
template< typename T >
inline void InitVec( Vec2<T>& inout, T arr[2] )
{
	inout(0)	= arr[0];
	inout(1)	= arr[1];
}


// Init Vector with zero
template< typename T >
void InitZero( Vec2<T>& inout )
{
	inout(0)	= 0;
	inout(1)	= 0;
}


// Reverse
template< typename T >
void Reverse( Vec2<T>& out, const Vec2<T>& in )
{
	out(0) = -in(0);
	out(1) = -in(1);
}


// Reverse
template< typename T >
void Reverse( Vec2<T>& inout )
{
	inout(0) = -inout(0);
	inout(1) = -inout(1);
}


// Add
template< typename T >
void Add( Vec2<T>& out, const Vec2<T>& in1, const Vec2<T>& in2 )
{
	out(0)	= in1(0) + in2(0);
	out(1)	= in1(1) + in2(1);
}


template< typename T >
inline void Add( Vec2<T>& inout, const Vec2<T>& in )
{
	inout(0) += in(0);
	inout(1) += in(1);
}


// Subtract
template< typename T >
void Subtract( Vec2<T> &out, const Vec2<T>& in1, const Vec2<T>& in2 )
{
	out(0)	= in1(0) - in2(0);
	out(1)	= in1(1) - in2(1);
}


template< typename T >
void Subtract( Vec2<T> &inout, const Vec2<T>& in )
{
	inout(0)	-= in(0);
	inout(1)	-= in(1);
}


// Dot product
template< typename T >
inline T DotProduct( const Vec2<T>& in1, const Vec2<T>& in2 )
{
	return in1(0) * in2(0) + in1(1) * in2(1);
}


// Cross product
template< typename T >
inline void CrossProduct( Vec2<T>& out, const Vec2<T>& in1, const Vec2<T>& in2 )
{
	out(0) = 0;
	out(1) = 0;
	out(2) = in1(0) * in2(1) - in1(1) * in2(0);
}


// Length
template< typename T >
inline T Length( const Vec2<T>& in )
{
	return sqrt( Max( in(0) * in(0) + in(1) * in(1), ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Length
template< typename T >
inline T LengthSqrd( const Vec2<T>& in )
{
	return in(0) * in(0) + in(1) * in(1);
}


// Distance between two vectors
template< typename T >
inline T Distance( const Vec2<T>& in1, const Vec2<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	return	sqrt( Max( dx * dx + dy * dy, ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Distance between two vectors
template< typename T >
inline T DistanceSqrd( const Vec2<T>& in1, const Vec2<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	return	dx * dx + dy * dy;
}


// Normalize
template< typename T >
inline void Normalize( Vec2<T>& inout )
{
	T length_inv	= ( T )1.0 / sqrt( Max( inout(0) * inout(0) + inout(1) * inout(1), ( std::numeric_limits<T>::min )( ) ) );
	inout(0) *= length_inv;
	inout(1) *= length_inv;
}


// Scale
template< typename T >
inline void Scale( Vec2<T>& inout, T scale )
{
	inout(0) *= scale;
	inout(1) *= scale;
}


template< typename T >
inline void Scale( Vec2<T>& out, const Vec2<T>& in, T scale )
{
	out(0) = in(0) * scale;
	out(1) = in(1) * scale;
}


template< typename T >
inline void Max( Vec2<T>& out, const Vec2<T>& in1, const Vec2<T>& in2 )
{
	out(0)	= in1(0) > in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) > in2(1) ? in1(1) : in2(1);
}


template< typename T >
inline void Min( Vec2<T>& out, const Vec2<T>& in1, const Vec2<T>& in2 )
{
	out(0)	= in1(0) < in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) < in2(1) ? in1(1) : in2(1);
}


// Clamp
template< typename T >
inline void Clamp( Vec2<T>& inout, const Vec2<T>& minVal, const Vec2<T>& maxVal )
{
	inout(0) = Max( Min( inout(0), maxVal(0) ), minVal(0) );
	inout(1) = Max( Min( inout(1), maxVal(1) ), minVal(1) );
}


// Lerp
template< typename T >
inline void Lerp( Vec2<T>& out, const Vec2<T>& start, const Vec2<T>& end, T percent )
{
	out(0)	= start(0) + percent * ( end(0) - start(0) );
	out(1)	= start(1) + percent * ( end(1) - start(1) );
}


// Spherilca Linear Interpolation
template< typename T >
inline void Slerp( Vec2<T>& out, const Vec2<T>& start, const Vec2<T>& end, T percent )
{
	// Dot product - the cosine of the angle between 2 vectors.
	T dot = (T)DotProduct( start, end );// Vector3.Dot(start, end);     
										// Clamp it to be in the range of Acos()
										// This may be unnecessary, but floating point
										// precision can be a fickle mistress.
	Clamp( dot, (T)-1, (T)1 ); //Mathf.Clamp(dot, -1.0f, 1.0f);
							   // Acos(dot) returns the angle between start and end,
							   // And multiplying that by percent returns the angle between
							   // start and the final result.
	T theta = (T)acos( dot ) * percent;//Mathf.Acos(dot)*percent;
	Vec2<T> RelativeVec; //Vector3 RelativeVec = end - start*dot;
	RelativeVec(0) = end(0) - dot * start(0);
	RelativeVec(1) = end(1) - dot * start(1);

	Normalize( RelativeVec );//RelativeVec.Normalize();     // Orthonormal basis
							 // The final result.
							 //return ((start*Mathf.Cos(theta)) + (RelativeVec*Mathf.Sin(theta)));
	T cos_theta = (T)cos( theta );
	T sin_theta = (T)sin( theta );
	out(0)	= start(0) * cos_theta +  RelativeVec(0) * sin_theta;
	out(1)	= start(1) * cos_theta +  RelativeVec(1) * sin_theta;
}


// Normalized Linear Interpolation
template< typename T >
inline void Nlerp( Vec2<T>& out, const Vec2<T>& start, const Vec2<T>& end, T percent )
{
	Lerp( out, start, end, percent );
	Normalize( out );
	//return Lerp( start, end, percent ).normalized();
}


template< typename T >
inline bool IsSame( const Vec2<T>& in1, const Vec2<T>& in2 )
{
	return in1(0)==in2(0) && in1(1)==in2(1);
}


template< typename T >
inline void AddScaled( Vec2<T>& out, float coeff1, const Vec2<T>& in1, float coeff2, const Vec2<T>& in2 )
{
	out(0)	= coeff1 * in1(0) + coeff2 * in2(0);
	out(1)	= coeff1 * in1(1) + coeff2 * in2(1);
}


template< typename T >
inline void AddScaled( Vec2<T>& out, const Vec2<T>& in1, float coeff2, const Vec2<T>& in2 )
{
	out(0)	= in1(0) + coeff2 * in2(0);
	out(1)	= in1(1) + coeff2 * in2(1);
}


template< typename T >
inline void AddScaled( Vec2<T>& inout, const Vec2<T>& in, const T scale )
{
	inout(0) += in(0) * scale;
	inout(1) += in(1) * scale;
}




//##############################################################################//
//									3D Vector									//
//##############################################################################//

//template< typename T >
//union Vec3
//{
//	T	xyz[3];
//	T	rgb[3];
//	struct { T x, y, z; };
//	struct { T r, g, b; };
//
//	Vec3()
//	{
//		x = 0;
//		y = 0;
//		z = 0;
//	}
//
//	Vec3( T x_, T y_, T z_ )
//	{
//		x = x_;
//		y = y_;
//		z = z_;
//	}
//
//
//	T& operator()( int i ) &
//	{
//		return xyz[i];
//	}
//
//	// Subscription operator. ( called by following cases: "T& a = Matrix<T>(10,10)(n)", "auto&& a = Matrix<T>(20,4)(n)" )
//	T operator()( int i ) const&&
//	{
//		return std::move(xyz[i]);
//	}
//
//};





// Init Vector
template< typename T >
void InitVec( Vec3<T>& inout, T x, T y, T z )
{
	inout(0)	= x;
	inout(1)	= y;
	inout(2)	= z;
}


// Init Vector
template< typename T >
inline void InitVec( Vec3<T>& inout, T arr[3] )
{
	inout(0)	= arr[0];
	inout(1)	= arr[1];
	inout(2)	= arr[2];
}


// Init Vector with zero
template< typename T >
void InitZero( Vec3<T>& inout )
{
	inout(0)	= 0;
	inout(1)	= 0;
	inout(2)	= 0;
}


// Reverse
template< typename T >
void Reverse( Vec3<T>& out, const Vec3<T>& in )
{
	out(0) = -in(0);
	out(1) = -in(1);
	out(2) = -in(2);
}


// Reverse
template< typename T >
void Reverse( Vec3<T>& inout )
{
	inout(0) = -inout(0);
	inout(1) = -inout(1);
	inout(2) = -inout(2);
}


// Add
template< typename T >
inline void Add( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0) = in1(0) + in2(0);
	out(1) = in1(1) + in2(1);
	out(2) = in1(2) + in2(2);
}


template< typename T >
inline void Add( Vec3<T>& inout, const Vec3<T>& in )
{
	inout(0) += in(0);
	inout(1) += in(1);
	inout(2) += in(2);
}


template< typename T >
inline void AddScaled( Vec3<T>& inout, const Vec3<T>& in, const T scale )
{
	inout(0) += in(0) * scale;
	inout(1) += in(1) * scale;
	inout(2) += in(2) * scale;
}


// Subtract
template< typename T >
inline void Subtract( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0) = in1(0) - in2(0);
	out(1) = in1(1) - in2(1);
	out(2) = in1(2) - in2(2);
}


// Multiply
template< typename T >
inline void Multiply( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0) = in1(0) * in2(0);
	out(1) = in1(1) * in2(1);
	out(2) = in1(2) * in2(2);
}


// Divide
template< typename T >
inline void Divide( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0) = in1(0) / in2(0);
	out(1) = in1(1) / in2(1);
	out(2) = in1(2) / in2(2);
}


// Dot product
template< typename T >
inline T DotProduct( const Vec3<T>& in1, const Vec3<T>& in2 )
{
	return in1(0) * in2(0) + in1(1) * in2(1) + in1(2) * in2(2);
}


// Cross product
template< typename T >
inline void CrossProduct( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0) = in1(1) * in2(2) - in1(2) * in2(1);
	out(1) = in1(2) * in2(0) - in1(0) * in2(2);
	out(2) = in1(0) * in2(1) - in1(1) * in2(0);
}


// Length
template< typename T >
inline T Length( const Vec3<T>& in )
{
	return sqrt( Max( in(0) * in(0) + in(1) * in(1) + in(2) * in(2), ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Length
template< typename T >
inline T LengthSqrd( const Vec3<T>& in )
{
	return in(0) * in(0) + in(1) * in(1) + in(2) * in(2);
}


// Distance between two vectors
template< typename T >
inline T Distance( const Vec3<T>& in1, const Vec3<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	const T dz	= in1(2) - in2(2);
	return	sqrt( Max( dx * dx + dy * dy + dz * dz, ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Distance between two vectors
template< typename T >
inline T DistanceSqrd( const Vec3<T>& in1, const Vec3<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	const T dz	= in1(2) - in2(2);
	return	dx * dx + dy * dy + dz * dz;
}


// Normalize
template< typename T >
inline void Normalize( Vec3<T>& inout )
{
	T length_inv	= ( T )1.0 / sqrt( Max( inout(0) * inout(0) + inout(1) * inout(1) + inout(2) * inout(2), ( std::numeric_limits<T>::min )( ) ) );
	inout(0) *= length_inv;
	inout(1) *= length_inv;
	inout(2) *= length_inv;
}


// Scale
template< typename T >
inline void Scale( Vec3<T>& inout, T scale )
{
	inout(0) *= scale;
	inout(1) *= scale;
	inout(2) *= scale;
}


template< typename T >
inline void Scale( Vec3<T>& out, const Vec3<T>& in, T scale )
{
	out(0) = in(0) * scale;
	out(1) = in(1) * scale;
	out(2) = in(2) * scale;
}


template< typename T >
inline void Max( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0)	= in1(0) > in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) > in2(1) ? in1(1) : in2(1);
	out(2)	= in1(2) > in2(2) ? in1(2) : in2(2);
}


template< typename T >
inline void Min( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 )
{
	out(0)	= in1(0) < in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) < in2(1) ? in1(1) : in2(1);
	out(2)	= in1(2) < in2(2) ? in1(2) : in2(2);
}


// Clamp
template< typename T >
inline void Clamp( Vec3<T>& inout, const Vec3<T>& minVal, const Vec3<T>& maxVal )
{
	inout(0) = Max( Min( inout(0), maxVal(0) ), minVal(0) );
	inout(1) = Max( Min( inout(1), maxVal(1) ), minVal(1) );
	inout(2) = Max( Min( inout(2), maxVal(2) ), minVal(2) );
}


// Lerp
template< typename T >
inline void Lerp( Vec3<T>& out, const Vec3<T>& start, const Vec3<T>& end, T percent )
{
	out(0)	= start(0) + percent * ( end(0) - start(0) );
	out(1)	= start(1) + percent * ( end(1) - start(1) );
	out(2)	= start(2) + percent * ( end(2) - start(2) );
}


// Spherilca Linear Interpolation
template< typename T >
inline void Slerp( Vec3<T>& out, const Vec3<T>& start, const Vec3<T>& end, T percent )
{
	// Dot product - the cosine of the angle between 2 vectors.
	T dot = (T)DotProduct( start, end );// Vector3.Dot(start, end);     
										// Clamp it to be in the range of Acos()
										// This may be unnecessary, but floating point
										// precision can be a fickle mistress.
	Clamp( dot, (T)-1, (T)1 ); //Mathf.Clamp(dot, -1.0f, 1.0f);
							   // Acos(dot) returns the angle between start and end,
							   // And multiplying that by percent returns the angle between
							   // start and the final result.
	T theta = (T)acos( dot ) * percent;//Mathf.Acos(dot)*percent;
	Vec3<T>	RelativeVec; //Vector3 RelativeVec = end - start*dot;
	RelativeVec(0) = end(0) - dot * start(0);
	RelativeVec(1) = end(1) - dot * start(1);
	RelativeVec(2) = end(2) - dot * start(2);

	Normalize( RelativeVec );//RelativeVec.Normalize();     // Orthonormal basis
							 // The final result.
							 //return ((start*Mathf.Cos(theta)) + (RelativeVec*Mathf.Sin(theta)));
	T cos_theta = (T)cos( theta );
	T sin_theta = (T)sin( theta );
	out(0)	= start(0) * cos_theta +  RelativeVec(0) * sin_theta;
	out(1)	= start(1) * cos_theta +  RelativeVec(1) * sin_theta;
	out(2)	= start(2) * cos_theta +  RelativeVec(2) * sin_theta;
}


// Normalized Linear Interpolation
template< typename T >
inline void Nlerp( Vec3<T>& out, const Vec3<T>& start, const Vec3<T>& end, T percent )
{
	Lerp( out, start, end, percent );
	Normalize( out );
	//return Lerp( start, end, percent ).normalized();
}


template< typename T >
inline bool IsSame( const Vec3<T>& in1, const Vec3<T>& in2 )
{
	return in1(0)==in2(0) && in1(1)==in2(1) && in1(2)==in2(2);
}


template< typename T >
inline void AddScaled( Vec3<T>& out, float coeff1, const Vec3<T>& in1, float coeff2, const Vec3<T>& in2 )
{
	out(0)	= coeff1 * in1(0) + coeff2 * in2(0);
	out(1)	= coeff1 * in1(1) + coeff2 * in2(1);
	out(2)	= coeff1 * in1(2) + coeff2 * in2(2);
}


template< typename T >
inline void AddScaled( Vec3<T>& out, const Vec3<T>& in1, float coeff2, const Vec3<T>& in2 )
{
	out(0)	= in1(0) + coeff2 * in2(0);
	out(1)	= in1(1) + coeff2 * in2(1);
	out(2)	= in1(2) + coeff2 * in2(2);
}



//##############################################################################//
//									4D Vector									//
//##############################################################################//

//template< typename T >
//union Vec4
//{
//	T	xyzw[4];
//	T	rgba[4];
//	struct { T x, y, z, w; };
//	struct { T r, g, b, a; };
//	struct { Vec3<T>xyz; T w; };
//
//	Vec4() : xyz()
//	{
//		//x = 0;
//		//y = 0;
//		//z = 0;
//		w = 0;
//	}
//
//	Vec4( T x_, T y_, T z_, T w_ )
//	{
//		x = x_;
//		y = y_;
//		z = z_;
//		w = w_;
//	}
//
//	T& operator()( int i ) &
//	{
//		return xyzw[i];
//	}
//
//	// Subscription operator. ( called by following cases: "T& a = Matrix<T>(10,10)(n)", "auto&& a = Matrix<T>(20,4)(n)" )
//	T operator()( int i ) const&&
//	{
//		return std::move(xyzw[i]);
//	}
//
//};




// Init Vector
template< typename T >
void InitVec( Vec4<T>& inout, T x, T y, T z, T w )
{
	inout(0)	= x;
	inout(1)	= y;
	inout(2)	= z;
	inout(3)	= w;
}


// Init Vector
template< typename T >
inline void InitVec( Vec4<T>& inout, T arr[4] )
{
	inout(0)	= arr[0];
	inout(1)	= arr[1];
	inout(2)	= arr[2];
	inout(3)	= arr[3];
}


// Init Vector with zero
template< typename T >
void InitZero( Vec4<T>& inout )
{
	inout(0)	= 0;
	inout(1)	= 0;
	inout(2)	= 0;
	inout(3)	= 0;
}


// Reverse
template< typename T >
void Reverse( Vec4<T>& out, const Vec4<T>& in )
{
	out(0) = -in(0);
	out(1) = -in(1);
	out(2) = -in(2);
	out(3) = -in(3);
}


// Reverse
template< typename T >
void Reverse( Vec4<T>& inout )
{
	inout(0) = -inout(0);
	inout(1) = -inout(1);
	inout(2) = -inout(2);
	inout(3) = -inout(3);
}


// Add
template< typename T >
inline void Add( Vec4<T>& out, const Vec4<T>& in1, const Vec4<T>& in2 )
{
	out(0) = in1(0) + in2(0);
	out(1) = in1(1) + in2(1);
	out(2) = in1(2) + in2(2);
	out(3) = in1(3) + in2(3);
}


template< typename T >
inline void Add( Vec4<T>& inout, const Vec4<T>& in )
{
	inout(0) += in(0);
	inout(1) += in(1);
	inout(2) += in(2);
	inout(3) += in(3);
}


// Subtract
template< typename T >
inline void Subtract( Vec4<T>& out, const Vec4<T>& in1, const Vec4<T>& in2 )
{
	out(0) = in1(0) - in2(0);
	out(1) = in1(1) - in2(1);
	out(2) = in1(2) - in2(2);
	out(3) = in1(3) - in2(3);
}


// Dot product
template< typename T >
inline T DotProduct( const Vec4<T>& in1, const Vec4<T>& in2 )
{
	return in1(0) * in2(0) + in1(1) * in2(1) + in1(2) * in2(2) + in1(3) * in2(3);
}


// Cross product
//template< typename T >
//inline void CrossProduct( Vec3<T>& out, const Vec3<T>& in1, const Vec3<T>& in2 ) 
//{
//	out(0) = in1(1) * in2(2) - in1(2) * in2(1);
//	out(1) = in1(2) * in2(0) - in1(0) * in2(2);
//	out(2) = in1(0) * in2(1) - in1(1) * in2(0);
//}


// Length
template< typename T >
inline T Length( const Vec4<T>& in )
{
	return sqrt( Max( in(0) * in(0) + in(1) * in(1) + in(2) * in(2) + in(3) * in(3), ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Length
template< typename T >
inline T LengthSqrd( const Vec4<T>& in )
{
	return in(0) * in(0) + in(1) * in(1) + in(2) * in(2) + in(3) * in(3);
}


// Distance between two vectors
template< typename T >
inline T Distance( const Vec4<T>& in1, const Vec4<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	const T dz	= in1(2) - in2(2);
	const T dw	= in1(3) - in2(3);
	return	sqrt( Max( dx * dx + dy * dy + dz * dz + dw * dw, ( std::numeric_limits<T>::min )( ) ) );
}


// Squared Distance between two vectors
template< typename T >
inline T DistanceSqrd( const Vec4<T>& in1, const Vec4<T>& in2 )
{
	const T dx	= in1(0) - in2(0);
	const T dy	= in1(1) - in2(1);
	const T dz	= in1(2) - in2(2);
	const T dw	= in1(3) - in2(3);
	return	dx * dx + dy * dy + dz * dz + dw* dw;
}


// Normalize
template< typename T >
inline void Normalize( Vec4<T>& inout )
{
	T length_inv	= ( T )1.0 / sqrt( Max( inout(0) * inout(0) + inout(1) * inout(1) + inout(2) * inout(2) + inout(3) * inout(3), ( std::numeric_limits<T>::min )( ) ) );
	inout(0) *= length_inv;
	inout(1) *= length_inv;
	inout(2) *= length_inv;
	inout(3) *= length_inv;
}


// Scale
template< typename T >
inline void Scale( Vec4<T>& inout, T scale )
{
	inout(0) *= scale;
	inout(1) *= scale;
	inout(2) *= scale;
	inout(3) *= scale;
}


template< typename T >
inline void Scale( Vec4<T>& out, const Vec4<T>& in, T scale )
{
	out(0) = in(0) * scale;
	out(1) = in(1) * scale;
	out(2) = in(2) * scale;
	out(3) = in(3) * scale;
}


template< typename T >
inline void Max( Vec4<T>& out, const Vec4<T>& in1, const Vec4<T>& in2 )
{
	out(0)	= in1(0) > in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) > in2(1) ? in1(1) : in2(1);
	out(2)	= in1(2) > in2(2) ? in1(2) : in2(2);
	out(3)	= in1(3) > in2(3) ? in1(3) : in2(3);
}


template< typename T >
inline void Min( Vec4<T>& out, const Vec4<T>& in1, const Vec4<T>& in2 )
{
	out(0)	= in1(0) < in2(0) ? in1(0) : in2(0);
	out(1)	= in1(1) < in2(1) ? in1(1) : in2(1);
	out(2)	= in1(2) < in2(2) ? in1(2) : in2(2);
	out(3)	= in1(3) < in2(3) ? in1(3) : in2(3);
}


// Clamp
template< typename T >
inline void Clamp( Vec4<T>& inout, const Vec4<T>& minVal, const Vec4<T>& maxVal )
{
	inout(0) = Max( Min( inout(0), maxVal(0) ), minVal(0) );
	inout(1) = Max( Min( inout(1), maxVal(1) ), minVal(1) );
	inout(2) = Max( Min( inout(2), maxVal(2) ), minVal(2) );
	inout(3) = Max( Min( inout(3), maxVal(3) ), minVal(3) );
}


// Lerp
template< typename T >
inline void Lerp( Vec4<T>& out, const Vec4<T>& start, const Vec4<T>& end, T percent )
{
	out(0)	= start(0) + percent * ( end(0) - start(0) );
	out(1)	= start(1) + percent * ( end(1) - start(1) );
	out(2)	= start(2) + percent * ( end(2) - start(2) );
	out(3)	= start(3) + percent * ( end(3) - start(3) );
}


// Spherilca Linear Interpolation
template< typename T >
inline void Slerp( Vec4<T>& out, const Vec4<T>& start, const Vec4<T>& end, T percent )
{
	// Dot product - the cosine of the angle between 2 vectors.
	T dot = (T)DotProduct( start, end );// Vector3.Dot(start, end);
										// Clamp it to be in the range of Acos()
										// This may be unnecessary, but floating point
										// precision can be a fickle mistress.
	Clamp( dot, (T)-1, (T)1 ); //Mathf.Clamp(dot, -1.0f, 1.0f);
							   // Acos(dot) returns the angle between start and end,
							   // And multiplying that by percent returns the angle between
							   // start and the final result.
	T theta = (T)acos( dot ) * percent;//Mathf.Acos(dot)*percent;
	Vec4<T>	RelativeVec; //Vector3 RelativeVec = end - start*dot;
	RelativeVec(0) = end(0) - dot * start(0);
	RelativeVec(1) = end(1) - dot * start(1);
	RelativeVec(2) = end(2) - dot * start(2);
	RelativeVec(3) = end(3) - dot * start(3);

	Normalize( RelativeVec );//RelativeVec.Normalize();     // Orthonormal basis
							 // The final result.
							 //return ((start*Mathf.Cos(theta)) + (RelativeVec*Mathf.Sin(theta)));
	T cos_theta = (T)cos( theta );
	T sin_theta = (T)sin( theta );
	out(0)	= start(0) * cos_theta +  RelativeVec(0) * sin_theta;
	out(1)	= start(1) * cos_theta +  RelativeVec(1) * sin_theta;
	out(2)	= start(2) * cos_theta +  RelativeVec(2) * sin_theta;
	out(3)	= start(3) * cos_theta +  RelativeVec(3) * sin_theta;
}


// Normalized Linear Interpolation
template< typename T >
inline void Nlerp( Vec4<T>& out, const Vec4<T>& start, const Vec4<T>& end, T percent )
{
	Lerp( out, start, end, percent );
	Normalize( out );
	//return Lerp( start, end, percent ).normalized();
}


template< typename T >
inline bool IsSame( const Vec4<T>& in1, const Vec4<T>& in2 )
{
	return in1(0)==in2(0) && in1(1)==in2(1) && in1(2)==in2(2) && in1(3)==in2(3);
}


template< typename T >
inline void AddScaled( Vec4<T>& out, float coeff1, const Vec4<T>& in1, float coeff2, const Vec4<T>& in2 )
{
	out(0)	= coeff1 * in1(0) + coeff2 * in2(0);
	out(1)	= coeff1 * in1(1) + coeff2 * in2(1);
	out(2)	= coeff1 * in1(2) + coeff2 * in2(2);
	out(3)	= coeff1 * in1(3) + coeff2 * in2(3);
}


template< typename T >
inline void AddScaled( Vec4<T>& out, const Vec4<T>& in1, float coeff2, const Vec4<T>& in2 )
{
	out(0)	= in1(0) + coeff2 * in2(0);
	out(1)	= in1(1) + coeff2 * in2(1);
	out(2)	= in1(2) + coeff2 * in2(2);
	out(3)	= in1(3) + coeff2 * in2(3);
}



//##############################################################################//
//									4x4 Matrix									//
//##############################################################################//

// 4x4 matrix union. left to right multiplication order
template< typename T>
union Mat4
{
	struct
	{
		T	m00, m01, m02, m03,
			m10, m11, m12, m13,
			m20, m21, m22, m23,
			m30, m31, m32, m33;
	};

	struct
	{
		Vec4<T>	mat[4];
	};

	//T		m[4][4];
	T		m[16];


	Mat4()
	{
		m00=0; m01=0; m02=0; m03=0;
		m10=0; m11=0; m12=0; m13=0;
		m20=0; m21=0; m22=0; m23=0;
		m30=0; m31=0; m32=0; m33=0;
	}

};




template< typename T >
inline void MatInit( Mat4<T>& inout,
	const T& m00, const T& m01, const T& m02, const T& m03,
	const T& m10, const T& m11, const T& m12, const T& m13,
	const T& m20, const T& m21, const T& m22, const T& m23,
	const T& m30, const T& m31, const T& m32, const T& m33 )
{
	inout.m00 = m00; inout.m01 = m01; inout.m02 = m02; inout.m03 = m03;
	inout.m10 = m10; inout.m11 = m11; inout.m12 = m12; inout.m13 = m13;
	inout.m20 = m20; inout.m21 = m21; inout.m22 = m22; inout.m23 = m23;
	inout.m30 = m30; inout.m31 = m31; inout.m32 = m32; inout.m33 = m33;
}



template< typename T >
inline void MatIdentity( Mat4<T>& inout )
{
	inout.m00 = 1; inout.m01 = 0; inout.m02 = 0; inout.m03 = 0;
	inout.m10 = 0; inout.m11 = 1; inout.m12 = 0; inout.m13 = 0;
	inout.m20 = 0; inout.m21 = 0; inout.m22 = 1; inout.m23 = 0;
	inout.m30 = 0; inout.m31 = 0; inout.m32 = 0; inout.m33 = 1;
}



template< typename T >
inline void MatZero( Mat4<T>& inout )
{
	inout.m00 = 0; inout.m01 = 0; inout.m02 = 0; inout.m03 = 0;
	inout.m10 = 0; inout.m11 = 0; inout.m12 = 0; inout.m13 = 0;
	inout.m20 = 0; inout.m21 = 0; inout.m22 = 0; inout.m23 = 0;
	inout.m30 = 0; inout.m31 = 0; inout.m32 = 0; inout.m33 = 0;
}



// inverse　
template< typename T >
inline void MatInverse( Mat4<T>& out, T& pDeterminant, const Mat4<T>& in )
{
	//pDeterminant =	(in.m00 * in.m11 * in.m22 * in.m33) + (in.m00 * in.m12 * in.m23 * in.m31) + (in.m00 * in.m13 * in.m21 * in.m32)
	//			+	(in.m01 * in.m10 * in.m23 * in.m32) + (in.m01 * in.m12 * in.m20 * in.m33) + (in.m01 * in.m13 * in.m22 * in.m30)
	//			+	(in.m02 * in.m10 * in.m21 * in.m33) + (in.m02 * in.m11 * in.m23 * in.m30) + (in.m02 * in.m13 * in.m20 * in.m31)
	//			+	(in.m03 * in.m10 * in.m22 * in.m31) + (in.m03 * in.m11 * in.m20 * in.m32) + (in.m03 * in.m12 * in.m21 * in.m30)
	//			-	(in.m00 * in.m11 * in.m23 * in.m32) - (in.m00 * in.m12 * in.m21 * in.m33) - (in.m00 * in.m13 * in.m22 * in.m31)
	//			-	(in.m01 * in.m10 * in.m22 * in.m33) - (in.m01 * in.m12 * in.m23 * in.m30) - (in.m01 * in.m13 * in.m20 * in.m32)
	//			-	(in.m02 * in.m10 * in.m23 * in.m31) - (in.m02 * in.m11 * in.m20 * in.m33) - (in.m02 * in.m13 * in.m21 * in.m30)
	//			-	(in.m03 * in.m10 * in.m21 * in.m32) - (in.m03 * in.m11 * in.m22 * in.m30) - (in.m03 * in.m12 * in.m20 * in.m31);

	//	
	//out.m00	= ( (in.m11 * in.m22 * in.m33) + (in.m12 * in.m23 * in.m31) + (in.m13 * in.m21 * in.m32) - (in.m11 * in.m23 * in.m32) - (in.m12 * in.m21 * in.m33) - (in.m13 * in.m22 * in.m31) ) / pDeterminant;
	//out.m01	= ( (in.m01 * in.m23 * in.m32) + (in.m02 * in.m21 * in.m33) + (in.m03 * in.m22 * in.m31) - (in.m01 * in.m22 * in.m33) - (in.m02 * in.m23 * in.m31) - (in.m03 * in.m21 * in.m32) ) / pDeterminant;
	//out.m02	= ( (in.m01 * in.m12 * in.m33) + (in.m02 * in.m13 * in.m31) + (in.m03 * in.m11 * in.m32) - (in.m01 * in.m13 * in.m32) - (in.m02 * in.m11 * in.m33) - (in.m03 * in.m12 * in.m31) ) / pDeterminant;
	//out.m03	= ( (in.m01 * in.m13 * in.m22) + (in.m02 * in.m11 * in.m23) + (in.m03 * in.m12 * in.m21) - (in.m01 * in.m12 * in.m23) - (in.m02 * in.m13 * in.m21) - (in.m03 * in.m11 * in.m22) ) / pDeterminant;
	//	
	//out.m10	= ( (in.m10 * in.m23 * in.m32) + (in.m12 * in.m20 * in.m33) + (in.m13 * in.m22 * in.m30) - (in.m10 * in.m22 * in.m33) - (in.m12 * in.m23 * in.m30) - (in.m13 * in.m20 * in.m32) ) / pDeterminant;
	//out.m11	= ( (in.m00 * in.m22 * in.m33) + (in.m02 * in.m23 * in.m30) + (in.m03 * in.m20 * in.m32) - (in.m00 * in.m23 * in.m32) - (in.m02 * in.m20 * in.m33) - (in.m03 * in.m22 * in.m30) ) / pDeterminant;
	//out.m12	= ( (in.m00 * in.m13 * in.m32) + (in.m02 * in.m10 * in.m33) + (in.m03 * in.m12 * in.m30) - (in.m00 * in.m12 * in.m33) - (in.m02 * in.m13 * in.m30) - (in.m03 * in.m10 * in.m32) ) / pDeterminant;
	//out.m13	= ( (in.m00 * in.m12 * in.m23) + (in.m02 * in.m13 * in.m20) + (in.m03 * in.m10 * in.m22) - (in.m00 * in.m13 * in.m22) - (in.m02 * in.m10 * in.m23) - (in.m03 * in.m12 * in.m20) ) / pDeterminant;
	//	
	//out.m20	= ( (in.m10 * in.m21 * in.m33) + (in.m11 * in.m23 * in.m30) + (in.m13 * in.m20 * in.m31) - (in.m10 * in.m23 * in.m31) - (in.m11 * in.m20 * in.m33) - (in.m13 * in.m21 * in.m30) ) / pDeterminant;
	//out.m21	= ( (in.m00 * in.m23 * in.m31) + (in.m01 * in.m20 * in.m33) + (in.m03 * in.m21 * in.m30) - (in.m00 * in.m21 * in.m33) - (in.m01 * in.m23 * in.m30) - (in.m03 * in.m20 * in.m31) ) / pDeterminant;
	//out.m22	= ( (in.m00 * in.m11 * in.m33) + (in.m01 * in.m13 * in.m30) + (in.m03 * in.m10 * in.m31) - (in.m00 * in.m13 * in.m31) - (in.m01 * in.m10 * in.m33) - (in.m03 * in.m11 * in.m30) ) / pDeterminant;
	//out.m23	= ( (in.m00 * in.m13 * in.m21) + (in.m01 * in.m10 * in.m23) + (in.m03 * in.m11 * in.m20) - (in.m00 * in.m11 * in.m23) - (in.m01 * in.m13 * in.m20) - (in.m03 * in.m10 * in.m21) ) / pDeterminant;

	//out.m30	= ( (in.m10 * in.m22 * in.m31) + (in.m11 * in.m20 * in.m32) + (in.m12 * in.m21 * in.m30) - (in.m10 * in.m21 * in.m32) - (in.m11 * in.m22 * in.m30) - (in.m12 * in.m20 * in.m31) ) / pDeterminant;
	//out.m31	= ( (in.m00 * in.m21 * in.m32) + (in.m01 * in.m22 * in.m30) + (in.m02 * in.m20 * in.m31) - (in.m00 * in.m22 * in.m31) - (in.m01 * in.m20 * in.m32) - (in.m02 * in.m21 * in.m30) ) / pDeterminant;
	//out.m32	= ( (in.m00 * in.m12 * in.m31) + (in.m01 * in.m10 * in.m32) + (in.m02 * in.m11 * in.m30) - (in.m00 * in.m11 * in.m32) - (in.m01 * in.m12 * in.m30) - (in.m02 * in.m10 * in.m31) ) / pDeterminant;
	//out.m33	= ( (in.m00 * in.m11 * in.m22) + (in.m01 * in.m12 * in.m20) + (in.m02 * in.m10 * in.m21) - (in.m00 * in.m12 * in.m21) - (in.m01 * in.m10 * in.m22) - (in.m02 * in.m11 * in.m20) ) / pDeterminant;
	//
	T s0	= in.m00 * in.m11 - in.m10 * in.m01;
	T s1	= in.m00 * in.m12 - in.m10 * in.m02;
	T s2	= in.m00 * in.m13 - in.m10 * in.m03;
	T s3	= in.m01 * in.m12 - in.m11 * in.m02;
	T s4	= in.m01 * in.m13 - in.m11 * in.m03;
	T s5	= in.m02 * in.m13 - in.m12 * in.m03;

	T c5	= in.m22 * in.m33 - in.m32 * in.m23;
	T c4	= in.m21 * in.m33 - in.m31 * in.m23;
	T c3	= in.m21 * in.m32 - in.m31 * in.m22;
	T c2	= in.m20 * in.m33 - in.m30 * in.m23;
	T c1	= in.m20 * in.m32 - in.m30 * in.m22;
	T c0	= in.m20 * in.m31 - in.m30 * in.m21;


	// Should check for 0 determinant
	pDeterminant	= ( s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0 );
	T invdet	= (T)1 / pDeterminant;

	out.m00	= ( in.m11 * c5 - in.m12 * c4 + in.m13 * c3 ) * invdet;
	out.m01	= ( -in.m01 * c5 + in.m02 * c4 - in.m03 * c3 ) * invdet;
	out.m02	= ( in.m31 * s5 - in.m32 * s4 + in.m33 * s3 ) * invdet;
	out.m03	= ( -in.m21 * s5 + in.m22 * s4 - in.m23 * s3 ) * invdet;

	out.m10	= ( -in.m10 * c5 + in.m12 * c2 - in.m13 * c1 ) * invdet;
	out.m11	= ( in.m00 * c5 - in.m02 * c2 + in.m03 * c1 ) * invdet;
	out.m12	= ( -in.m30 * s5 + in.m32 * s2 - in.m33 * s1 ) * invdet;
	out.m13	= ( in.m20 * s5 - in.m22 * s2 + in.m23 * s1 ) * invdet;

	out.m20	= ( in.m10 * c4 - in.m11 * c2 + in.m13 * c0 ) * invdet;
	out.m21	= ( -in.m00 * c4 + in.m01 * c2 - in.m03 * c0 ) * invdet;
	out.m22	= ( in.m30 * s4 - in.m31 * s2 + in.m33 * s0 ) * invdet;
	out.m23	= ( -in.m20 * s4 + in.m21 * s2 - in.m23 * s0 ) * invdet;

	out.m30	= ( -in.m10 * c3 + in.m11 * c1 - in.m12 * c0 ) * invdet;
	out.m31	= ( in.m00 * c3 - in.m01 * c1 + in.m02 * c0 ) * invdet;
	out.m32	= ( -in.m30 * s3 + in.m31 * s1 - in.m32 * s0 ) * invdet;
	out.m33	= ( in.m20 * s3 - in.m21 * s1 + in.m22 * s0 ) * invdet;


}



template< typename T >
inline void MatInverse( Mat4<T>& out, const Mat4<T>& in )
{
	T s0	= in.m00 * in.m11 - in.m10 * in.m01;
	T s1	= in.m00 * in.m12 - in.m10 * in.m02;
	T s2	= in.m00 * in.m13 - in.m10 * in.m03;
	T s3	= in.m01 * in.m12 - in.m11 * in.m02;
	T s4	= in.m01 * in.m13 - in.m11 * in.m03;
	T s5	= in.m02 * in.m13 - in.m12 * in.m03;

	T c5	= in.m22 * in.m33 - in.m32 * in.m23;
	T c4	= in.m21 * in.m33 - in.m31 * in.m23;
	T c3	= in.m21 * in.m32 - in.m31 * in.m22;
	T c2	= in.m20 * in.m33 - in.m30 * in.m23;
	T c1	= in.m20 * in.m32 - in.m30 * in.m22;
	T c0	= in.m20 * in.m31 - in.m30 * in.m21;


	// Should check for 0 determinant
	T invdet	= (T)1 / ( s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0 );

	out.m00	= ( in.m11 * c5 - in.m12 * c4 + in.m13 * c3 ) * invdet;
	out.m01	= ( -in.m01 * c5 + in.m02 * c4 - in.m03 * c3 ) * invdet;
	out.m02	= ( in.m31 * s5 - in.m32 * s4 + in.m33 * s3 ) * invdet;
	out.m03	= ( -in.m21 * s5 + in.m22 * s4 - in.m23 * s3 ) * invdet;

	out.m10	= ( -in.m10 * c5 + in.m12 * c2 - in.m13 * c1 ) * invdet;
	out.m11	= ( in.m00 * c5 - in.m02 * c2 + in.m03 * c1 ) * invdet;
	out.m12	= ( -in.m30 * s5 + in.m32 * s2 - in.m33 * s1 ) * invdet;
	out.m13	= ( in.m20 * s5 - in.m22 * s2 + in.m23 * s1 ) * invdet;

	out.m20	= ( in.m10 * c4 - in.m11 * c2 + in.m13 * c0 ) * invdet;
	out.m21	= ( -in.m00 * c4 + in.m01 * c2 - in.m03 * c0 ) * invdet;
	out.m22	= ( in.m30 * s4 - in.m31 * s2 + in.m33 * s0 ) * invdet;
	out.m23	= ( -in.m20 * s4 + in.m21 * s2 - in.m23 * s0 ) * invdet;

	out.m30	= ( -in.m10 * c3 + in.m11 * c1 - in.m12 * c0 ) * invdet;
	out.m31	= ( in.m00 * c3 - in.m01 * c1 + in.m02 * c0 ) * invdet;
	out.m32	= ( -in.m30 * s3 + in.m31 * s1 - in.m32 * s0 ) * invdet;
	out.m33	= ( in.m20 * s3 - in.m21 * s1 + in.m22 * s0 ) * invdet;
}






// transpose
template< typename T >
inline void MatTranspose( Mat4<T>& out, const Mat4<T>& in )
{
	out.m00	= in.m00;
	out.m01	= in.m10;
	out.m02	= in.m20;
	out.m03	= in.m30;

	out.m10	= in.m01;
	out.m11	= in.m11;
	out.m12	= in.m21;
	out.m13	= in.m31;

	out.m20	= in.m02;
	out.m21	= in.m12;
	out.m22	= in.m22;
	out.m23	= in.m32;

	out.m30	= in.m03;
	out.m31	= in.m13;
	out.m32	= in.m23;
	out.m33	= in.m33;
}


// add
template< typename T >
inline void Add( Mat4<T>& out, const Mat4<T>& in1, const Mat4<T>& in2 )
{
	out.m00	= in1.m00 + in2.m00;
	out.m01	= in1.m01 + in2.m01;
	out.m02	= in1.m02 + in2.m02;
	out.m03	= in1.m03 + in2.m03;

	out.m10	= in1.m10 + in2.m10;
	out.m11	= in1.m11 + in2.m11;
	out.m12	= in1.m12 + in2.m12;
	out.m13	= in1.m13 + in2.m13;

	out.m20	= in1.m20 + in2.m20;
	out.m21	= in1.m21 + in2.m21;
	out.m22	= in1.m22 + in2.m22;
	out.m23	= in1.m23 + in2.m23;

	out.m30	= in1.m30 + in2.m30;
	out.m31	= in1.m31 + in2.m31;
	out.m32	= in1.m32 + in2.m32;
	out.m33	= in1.m33 + in2.m33;
}


// subtract
template< typename T >
inline void Subtract( Mat4<T>& out, const Mat4<T>& in1, const Mat4<T>& in2 )
{
	out.m00	= in1.m00 - in2.m00;
	out.m01	= in1.m01 - in2.m01;
	out.m02	= in1.m02 - in2.m02;
	out.m03	= in1.m03 - in2.m03;

	out.m10	= in1.m10 - in2.m10;
	out.m11	= in1.m11 - in2.m11;
	out.m12	= in1.m12 - in2.m12;
	out.m13	= in1.m13 - in2.m13;

	out.m20	= in1.m20 - in2.m20;
	out.m21	= in1.m21 - in2.m21;
	out.m22	= in1.m22 - in2.m22;
	out.m23	= in1.m23 - in2.m23;

	out.m30	= in1.m30 - in2.m30;
	out.m31	= in1.m31 - in2.m31;
	out.m32	= in1.m32 - in2.m32;
	out.m33	= in1.m33 - in2.m33;
}


// note: in(3) is assumed to be 1.0
template< typename T >
inline void Multiply( Vec3<T>& out, const Mat4<T>& mat, const Vec3<T>& in )
{
	out(0)	= mat.m00 * in(0) +  mat.m01 * in(1) +  mat.m02 * in(2) + mat.m03;
	out(1)	= mat.m10 * in(0) +  mat.m11 * in(1) +  mat.m12 * in(2) + mat.m13;
	out(2)	= mat.m20 * in(0) +  mat.m21 * in(1) +  mat.m22 * in(2) + mat.m23;
}


template< typename T >
inline void Multiply( Vec4<T>& out, const Mat4<T>& mat, const Vec3<T>& in )
{
	out(0)	= mat.m00 * in(0) +  mat.m01 * in(1) +  mat.m02 * in(2) + mat.m03;
	out(1)	= mat.m10 * in(0) +  mat.m11 * in(1) +  mat.m12 * in(2) + mat.m13;
	out(2)	= mat.m20 * in(0) +  mat.m21 * in(1) +  mat.m22 * in(2) + mat.m23;
	out(3)	= mat.m30 * in(0) +  mat.m31 * in(1) +  mat.m32 * in(2) + mat.m33;
}



template< typename T >
inline void Multiply( Vec4<T>& out, const Mat4<T>& mat, const Vec4<T>& in )
{
	out(0)	= mat.m00 * in(0) +  mat.m01 * in(1) +  mat.m02 * in(2) + mat.m03 * in(3);
	out(1)	= mat.m10 * in(0) +  mat.m11 * in(1) +  mat.m12 * in(2) + mat.m13 * in(3);
	out(2)	= mat.m20 * in(0) +  mat.m21 * in(1) +  mat.m22 * in(2) + mat.m23 * in(3);
	out(3)	= mat.m30 * in(0) +  mat.m31 * in(1) +  mat.m32 * in(2) + mat.m33 * in(3);
}


template< typename T >
inline void Multiply( Mat4<T>& out, const Mat4<T>& in1, const Mat4<T>& in2 )
{
	out.m00 = in1.m00*in2.m00 + in1.m01*in2.m10 + in1.m02*in2.m20 + in1.m03*in2.m30;
	out.m01 = in1.m00*in2.m01 + in1.m01*in2.m11 + in1.m02*in2.m21 + in1.m03*in2.m31;
	out.m02 = in1.m00*in2.m02 + in1.m01*in2.m12 + in1.m02*in2.m22 + in1.m03*in2.m32;
	out.m03 = in1.m00*in2.m03 + in1.m01*in2.m13 + in1.m02*in2.m23 + in1.m03*in2.m33;

	out.m10 = in1.m10*in2.m00 + in1.m11*in2.m10 + in1.m12*in2.m20 + in1.m13*in2.m30;
	out.m11 = in1.m10*in2.m01 + in1.m11*in2.m11 + in1.m12*in2.m21 + in1.m13*in2.m31;
	out.m12 = in1.m10*in2.m02 + in1.m11*in2.m12 + in1.m12*in2.m22 + in1.m13*in2.m32;
	out.m13 = in1.m10*in2.m03 + in1.m11*in2.m13 + in1.m12*in2.m23 + in1.m13*in2.m33;

	out.m20 = in1.m20*in2.m00 + in1.m21*in2.m10 + in1.m22*in2.m20 + in1.m23*in2.m30;
	out.m21 = in1.m20*in2.m01 + in1.m21*in2.m11 + in1.m22*in2.m21 + in1.m23*in2.m31;
	out.m22 = in1.m20*in2.m02 + in1.m21*in2.m12 + in1.m22*in2.m22 + in1.m23*in2.m32;
	out.m23 = in1.m20*in2.m03 + in1.m21*in2.m13 + in1.m22*in2.m23 + in1.m23*in2.m33;

	out.m30 = in1.m30*in2.m00 + in1.m31*in2.m10 + in1.m32*in2.m20 + in1.m33*in2.m30;
	out.m31 = in1.m30*in2.m01 + in1.m31*in2.m11 + in1.m32*in2.m21 + in1.m33*in2.m31;
	out.m32 = in1.m30*in2.m02 + in1.m31*in2.m12 + in1.m32*in2.m22 + in1.m33*in2.m32;
	out.m33 = in1.m30*in2.m03 + in1.m31*in2.m13 + in1.m32*in2.m23 + in1.m33*in2.m33;
}






template< typename T >
inline void MatScale( Mat4<T>& out, const Vec3<T>& scale )
{
	out.m00 = scale(0);	out.m01 = 0;		out.m02 = 0;		out.m03 = 0;
	out.m10 = 0;		out.m11 = scale(1);	out.m12 = 0;		out.m13 = 0;
	out.m20 = 0;		out.m21 = 0;		out.m22 = scale(2);	out.m23 = 0;
	out.m30 = 0;		out.m31 = 0;		out.m32 = 0;		out.m33 = 1;
}


template< typename T >
inline void MatScale( Mat4<T>& out, T sx, T sy, T sz )
{
	out.m00 = sx;		out.m01 = 0;		out.m02 = 0;		out.m03 = 0;
	out.m10 = 0;		out.m11 = sy;		out.m12 = 0;		out.m13 = 0;
	out.m20 = 0;		out.m21 = 0;		out.m22 = sz;		out.m23 = 0;
	out.m30 = 0;		out.m31 = 0;		out.m32 = 0;		out.m33 = 1;
}


template< typename T >
inline void MatTranslation( Mat4<T>& out, const Vec3<T>& vec )
{
	out.m00 = 1; out.m01 = 0; out.m02 = 0; out.m03 = vec(0);
	out.m10 = 0; out.m11 = 1; out.m12 = 0; out.m13 = vec(1);
	out.m20 = 0; out.m21 = 0; out.m22 = 1; out.m23 = vec(2);
	out.m30 = 0; out.m31 = 0; out.m32 = 0; out.m33 = 1;
}


template< typename T >
inline void MatTranslation( Mat4<T>& out, T vx, T vy, T vz )
{
	out.m00 = 1; out.m01 = 0; out.m02 = 0; out.m03 = vx;
	out.m10 = 0; out.m11 = 1; out.m12 = 0; out.m13 = vy;
	out.m20 = 0; out.m21 = 0; out.m22 = 1; out.m23 = vz;
	out.m30 = 0; out.m31 = 0; out.m32 = 0; out.m33 = 1;
}


template< typename T >
inline void MatRotationX( Mat4<T>& mat, T theta )
{
	mat.m00 = 1;		mat.m01 = 0;			mat.m02 = 0;			mat.m03 = 0;
	mat.m10 = 0;		mat.m11 = cos( theta );	mat.m12 = -sin( theta );	mat.m13 = 0;
	mat.m20 = 0;		mat.m21 = sin( theta );	mat.m22 = cos( theta );	mat.m23 = 0;
	mat.m30 = 0;		mat.m31 = 0;			mat.m32 = 0;			mat.m33 = 1;
}


template< typename T >
inline void MatRotationY( Mat4<T>& mat, T theta )
{
	mat.m00 = cos( theta );	mat.m01 = 0;		mat.m02 = sin( theta );	mat.m03 = 0;
	mat.m10 = 0;			mat.m11 = 1;		mat.m12 = 0;			mat.m13 = 0;
	mat.m20 = -sin( theta );	mat.m21 = 0;		mat.m22 = cos( theta );	mat.m23 = 0;
	mat.m30 = 0;			mat.m31 = 0;		mat.m32 = 0;			mat.m33 = 1;
}


template< typename T >
inline void MatRotationZ( Mat4<T>& mat, T theta )
{
	mat.m00 = cos( theta );	mat.m01 = -sin( theta );	mat.m02 = 0;		mat.m03 = 0;
	mat.m10 = sin( theta );	mat.m11 = cos( theta );	mat.m12 = 0;		mat.m13 = 0;
	mat.m20 = 0;			mat.m21 = 0;			mat.m22 = 1;		mat.m23 = 0;
	mat.m30 = 0;			mat.m31 = 0;			mat.m32 = 0;		mat.m33 = 1;
}


// 座標変換行列を作成する(右手座標系/左手座標系の違いに注意！).行順
// u: horizontal, v: vertical, n: forward
template< typename T >
inline void MatViewGL( Mat4<T>& out, const Vec3<T>& u, const Vec3<T>& v, const Vec3<T>& n, const Vec3<T>& c )
{
	out.m00 = -u(0);		out.m01 = -u(1);		out.m02 = -u(2);		out.m03 = DotProduct( u, c );
	out.m10 = v(0);		out.m11 = v(1);		out.m12 = v(2);		out.m13 = -DotProduct( v, c );
	out.m20 = -n(0);		out.m21 = -n(1);		out.m22 = -n(2);		out.m23 = DotProduct( n, c );
	out.m30 = 0;		out.m31 = 0;		out.m32 = 0;		out.m33 = 1;
}


template< typename T >
inline void MatView( Mat4<T>& out, const Vec3<T>& u, const Vec3<T>& v, const Vec3<T>& n, const Vec3<T>& pos )
{
	out.m00 = u(0);	out.m01 = u(1);	out.m02 = u(2);	out.m03 = DotProduct( u, pos );
	out.m10 = v(0);	out.m11 = v(1);	out.m12 = v(2);	out.m13 = DotProduct( v, pos );
	out.m20 = n(0);	out.m21 = n(1);	out.m22 = n(2);	out.m23 = DotProduct( n, pos );
	out.m30 = 0;	out.m31 = 0;	out.m32 = 0;	out.m33 = 1;
}



// 射影変換行列を作成する(gluPerspectiveと同じ).行順
template< typename T >
inline void MatPerspectiveFov( Mat4<T>& out, T fovy, T aspect, T znear, T zfar )
{
	T	depth = (znear)-( zfar );
	T	f = ( T )1.0 / tan( ( T )0.5*( fovy ) );

	out.m00 = f / ( aspect );		out.m01 = 0;			out.m02 = 0;							out.m03 = 0;
	out.m10 = 0;				out.m11 = f;			out.m12 = 0;							out.m13 = 0;
	out.m20 = 0;				out.m21 = 0;			out.m22 = ( (zfar)+( znear ) ) / depth;		out.m23 = (T)2*( zfar )*( znear ) / depth;
	out.m30 = 0;				out.m31 = 0;			out.m32 = -1;							out.m33 = 0;
}


// 射影変換行列を作成する(glFrustumと同じ).行順
template< typename T >
inline void MatPerspectiveOffCenter( Mat4<T>& out, T left, T right, T bottom, T top, T znear, T zfar )
{
	out.m00 = 2*znear/( right-left );	out.m01 = 0;						out.m02 = ( right+left )/( right-left );	out.m03 = 0;
	out.m10 = 0;					out.m11 = 2*znear/( top-bottom );		out.m12 = ( top+bottom )/( top-bottom );	out.m13 = 0;
	out.m20 = 0;					out.m21 = 0;						out.m22 = -( zfar+znear )/( zfar-znear );	out.m23 = -2*zfar*znear/( zfar-znear ); //znear*zfar/(znear-zfar);
	out.m30 = 0;					out.m31 = 0;						out.m32 = -1;							out.m33 = 0;
}



// 射影変換行列を作成する(glOrthoと等価).
template< typename T >
inline void MatOrtho( Mat4<T>& out, T left, T right, T bottom, T top, T znear, T zfar )
{
	T width		= right -left;
	T height	= top - bottom;
	T depth		= zfar - znear;

	out.m00 = (T)2 / width;		out.m01 = 0;				out.m02 = 0;				out.m03 = -( right + left ) / width;
	out.m10 = 0;				out.m11 = (T)2 / height;	out.m12 = 0;				out.m13 = -( top + bottom ) / height;
	out.m20 = 0;				out.m21 = 0;				out.m22 = (T)-2 / depth;	out.m23 = -( zfar + znear ) / depth;
	out.m30 = 0;				out.m31 = 0;				out.m32 = 0;				out.m33 = 1;
}




//##############################################################################//
//									4x4 Matrix									//
//##############################################################################//







inline void InitQuat( Vec4f& quat, float angle, float x, float y, float z )
{
	const float theta = angle * 0.5f;
	const float sine_theta = sin( theta );
	quat(3)	= cos( theta );
	quat(0)	= x * sine_theta;
	quat(1)	= y * sine_theta;
	quat(2)	= z * sine_theta;
}


inline float length( Vec4f quat )
{
	return sqrt( quat(0)*quat(0) + quat(1)*quat(1) +quat(2)*quat(2) + quat(3)*quat(3) );
}


inline Vec4f normalize( Vec4f quat )
{
	float L_inv = 1.0f / length( quat );

	quat(0) *= L_inv;
	quat(1) *= L_inv;
	quat(2) *= L_inv;
	quat(3) *= L_inv;

	return quat;
}


inline Vec4f conjugate( Vec4f quat )
{
	quat(0) = -quat(0);
	quat(1) = -quat(1);
	quat(2) = -quat(2);
	return quat;
}


// クォータニオン合成蓄積誤差除去
inline void rm_cal_err_qt( Vec4f quat )
{

	float s = length( quat );
	if( s <= 0.0f ) return;
	quat(0) /= s;
	quat(1) /= s;
	quat(2) /= s;
	quat(3) /= s;
}



inline Vec4f mult( Vec4f A, Vec4f B )
{
	Vec4f C;

	C(0) = A(3)*B(0) + A(0)*B(3) + A(1)*B(2) - A(2)*B(1);
	C(1) = A(3)*B(1) - A(0)*B(2) + A(1)*B(3) + A(2)*B(0);
	C(2) = A(3)*B(2) + A(0)*B(1) - A(1)*B(0) + A(2)*B(3);
	C(3) = A(3)*B(3) - A(0)*B(0) - A(1)*B(1) - A(2)*B(2);


	return C;
}


// *=も大丈夫なバージョン
//inline Vec4f mult( Vec4f lpP, Vec4f lpQ )
//{
//	Vec4f lpR;
//
//    float pw, px, py, pz;
//    float qw, qx, qy, qz;
//
//    pw = lpP(3); px = lpP(0); py = lpP(1); pz = lpP(2);
//    qw = lpQ(3); qx = lpQ(0); qy = lpQ(1); qz = lpQ(2);
//
//    lpR(3) = pw * qw - px * qx - py * qy - pz * qz;
//    lpR(0) = pw * qx + px * qw + py * qz - pz * qy;
//    lpR(1) = pw * qy - px * qz + py * qw + pz * qx;
//    lpR(2) = pw * qz + px * qy - py * qx + pz * qw;
//
//	return lpR;
//}












//template< typename T >
//union Quaternion
//{
//	struct { T w, x, y, z; };
//};



template< typename T >
inline void InitQuat( Quaternion<T>& quat, T angle, T x, T y, T z )
{
	const double theta	= (double)angle * 0.5;
	const T sine_theta	= (T)sin( theta );

	Vec3f axis			={ x, y, z };
	Normalize( axis );

	quat(3)	= (T)cos( theta );
	quat(0)	= axis(0) * sine_theta;
	quat(1)	= axis(1) * sine_theta;
	quat(2)	= axis(2) * sine_theta;
}


template< typename T >
inline void InitQuat( Quaternion<T>& quat, T angle, const Vec3<T>& axis )
{
	const double theta = (double)angle * 0.5;
	const T sine_theta = (T)sin( theta );
	quat(3)	= (T)cos( theta );
	quat(0)	= axis(0) * sine_theta;
	quat(1)	= axis(1) * sine_theta;
	quat(2)	= axis(2) * sine_theta;
}


template< typename T >
inline void QuatIdentity( Quaternion<T>& quat )
{
	quat(3)	= (T)1;
	quat(0)	= 0;
	quat(1)	= 0;
	quat(2)	= 0;
}


template< typename T >
inline void QuatConjugate( Quaternion<T>& quat )
{
	quat(0) = -quat(0);
	quat(1) = -quat(1);
	quat(2) = -quat(2);
}


template< typename T >
inline void QuatConjugate( Quaternion<T>& out, const Quaternion<T>& quat )
{
	out(3)	= quat(3);
	out(0)	= -quat(0);
	out(1)	= -quat(1);
	out(2)	= -quat(2);
}


template< typename T >
inline T Length( const Quaternion<T>& quat )
{
	return sqrt( quat(3)*quat(3) + quat(0)*quat(0) + quat(1)*quat(1) +quat(2)*quat(2) );
}


template < typename T >
inline void Normalize( Quaternion<T>& quat )
{
	float L_inv = ( T )1.0 / Length( quat );

	quat(3) *= L_inv;
	quat(0) *= L_inv;
	quat(1) *= L_inv;
	quat(2) *= L_inv;
}


// クォータニオン合成蓄積誤差除去
template< typename T >
inline void rm_cal_err_qt( Quaternion<T>& quat )
{
	T s = Length( quat );
	if( s <= (T)0 ) return;

	quat(3) /= s;
	quat(0) /= s;
	quat(1) /= s;
	quat(2) /= s;
}


// Quaternionの乗算. Aの右にBをかける.B回転が最初に作用して、次にA回転が作用する
template< typename T >
inline void Multiply( Quaternion<T>& C, const Quaternion<T>& A, const Quaternion<T>& B )
{
	C(0) = A(3)*B(0) + A(0)*B(3) + A(1)*B(2) - A(2)*B(1);
	C(1) = A(3)*B(1) - A(0)*B(2) + A(1)*B(3) + A(2)*B(0);
	C(2) = A(3)*B(2) + A(0)*B(1) - A(1)*B(0) + A(2)*B(3);
	C(3) = A(3)*B(3) - A(0)*B(0) - A(1)*B(1) - A(2)*B(2);
}

// Quaternionの乗算. lpPの右にlpQをかける.lpQ回転が最初に作用して、次にlpP回転が作用する
template< typename T >
inline void Multiply_( Quaternion<T>& lpR, const Quaternion<T>& lpP, const Quaternion<T>& lpQ )
{
	T pw, px, py, pz;
	T qw, qx, qy, qz;

	pw = lpP(3); px = lpP(0); py = lpP(1); pz = lpP(2);
	qw = lpQ(3); qx = lpQ(0); qy = lpQ(1); qz = lpQ(2);

	lpR(3) = pw * qw - px * qx - py * qy - pz * qz;
	lpR(0) = pw * qx + px * qw + py * qz - pz * qy;
	lpR(1) = pw * qy - px * qz + py * qw + pz * qx;
	lpR(2) = pw * qz + px * qy - py * qx + pz * qw;
}



template< typename T >
inline void Rotate( Vec3<T>& inout, const Quaternion<T>& quat )
{
	Quaternion<T>	q1, quatConjugate, result;
	Quaternion<T>	quat_inout ={ 0, inout(0), inout(1), inout(2) };

	QuatConjugate( quatConjugate, quat );

	// quatConjugate * quat_inout * quat( TODO: 注意!. 最後にかけたものが最初に有効になる。quatConjugate * inout * quatは後ろからかける)
	Multiply( q1, quat, quat_inout );
	Multiply( result, q1, quatConjugate );

	inout(0) = result(0);
	inout(1) = result(1);
	inout(2) = result(2);
}


// 姿勢変化をクォータニオンに変換する.Z軸方向がforwardに、Y軸方向がupにアライメントされる
template< typename T >
inline void QuaternionLookAt( Quaternion<T>& out, const Vec3<T>& forward, const Vec3<T>& up, const Vec3<T>& basis_forward, const Vec3<T>& basis_up )
{
	//const Vec3<T> Y_AXIS = { 0.0, 1.0, 0.0 };
	//const Vec3<T> Z_AXIS = { 0.0, 0.0, 1.0 };

	//=============== Z軸からforwardへの回転のクォータニオンを計算する ================//
	// Z_AXISとforwardの外積を計算して、回転軸を取得する
	Vec3<T> rot_axis_forward;
	CrossProduct( rot_axis_forward, forward, basis_forward/*Z_AXIS*/ );
	Normalize( rot_axis_forward );

	// Z_AXISとforwardの内積を計算して、回転角を取得する( TODO: 注意!. 外積の回転方向は右ねじの法則。クォータニオンは逆)
	T angle_forward	= -acos( DotProduct( basis_forward/*Z_AXIS*/, forward ) );

	// forwardのクォータニオンを作成する
	Quaternion<T> quat_forward;
	InitQuat( quat_forward, angle_forward, rot_axis_forward );


	//=============== forward分のクォータニオンをY_AXIS軸に作用させる =================//
	Vec3<T> up_before = basis_up;//Y_AXIS;
	Rotate( up_before, quat_forward );

	//================ up_beforeからupへのクォータニオンを計算する ===================//
	// 回転軸を計算する
	Vec3<T> rot_axis_up;
	CrossProduct( rot_axis_up, up, up_before );
	Normalize( rot_axis_up );

	// up_beforeとupの内積を計算して、回転角を取得する( TODO: 注意!. 外積の回転方向は右ねじの法則。クォータニオンは逆)
	T angle_up	= -acos( DotProduct( up_before, up ) );

	// upのクォータニオンを作成する
	Quaternion<T> quat_up;
	InitQuat( quat_up, angle_up, rot_axis_up );


	//========================== クォータニオンを合成する ============================//
	Multiply( out, quat_up, quat_forward );
	//out = quat_forward;
}




// クォータニオンから回転行列への変換
template< typename T >
inline void Quat2Mat( Mat4<T>& out, const Quaternion<T>& quat )
{
	T x2_2 = quat(0) * quat(0) * 2;
	T xy_2 = quat(0) * quat(1) * 2;
	T xz_2 = quat(0) * quat(2) * 2;
	T wx_2 = quat(0) * quat(3) * 2;

	T y2_2 = quat(1) * quat(1) * 2;
	T yz_2 = quat(1) * quat(2) * 2;
	T wy_2 = quat(1) * quat(3) * 2;

	T z2_2 = quat(2) * quat(2) * 2;
	T zw_2 = quat(2) * quat(3) * 2;

	out.m00	= 1-y2_2-z2_2;	out.m01 = xy_2-zw_2;	out.m02 = xz_2+wy_2;	out.m03 = 0;
	out.m10 = xy_2 + zw_2;	out.m11 = 1-x2_2-z2_2;	out.m12 = yz_2-wx_2;	out.m13 = 0;
	out.m20 = xz_2 - wy_2;	out.m21 = yz_2 + wx_2;	out.m22 = 1-x2_2-y2_2;	out.m23 = 0;
	out.m30 = 0;			out.m31 = 0;			out.m32 = 0;			out.m33 = 1;
}




#endif // VECTOR_OPERATIONS_H