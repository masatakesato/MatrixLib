#ifndef VECTOR_OPERATIONS_H
#define	VECTOR_OPERATIONS_H


#include	<math.h>



// Euclidean Distance
template< typename T>
T EuclideanDistance( int dim, const T vec1[], const T vec2[] )
{
	double dist = 0;
	for( int i=0; i<dim; ++i )
	{
		dist += double( ( vec1[i] - vec2[i] ) * ( vec1[i] - vec2[i] ) );
	}
	return (T)sqrt( dist );
}


// Chebyshev Distance
template< typename T>
T ChebyshevDistance( int dim, const T vec1[], const T vec2[] )
{
	double dist = 0;
	for( int i=0; i<dim; ++i )
	{
		T d = fabs( double( vec1[i] - vec2[i] ) );
		dist = d > dist ? d : dist;
	}
	return dist;
}


// Standardized Euclidean distance



// Manhattan Distance
template< typename T>
T MahnattannDistance( int dim, const T vec1[], const T vec2[] )
{
	double dist = 0;
	for( int i=0; i<dim; ++i )
	{
		dist += fabs( double( vec1[i] - vec2[i] ) );
	}
	return (T)dist;
}





// Cosine similarity
template< typename T>
static T CosineSimilarity( int dim, const T vec1[], const T vec2[] )
{
	double inner_product = 0;
	double len1	= 0;
	double len2 = 0;

	for( int i=0; i<dim; ++i )
	{
		inner_product += double( vec1[i] * vec2[i] );
		len1 += double( vec1[i] * vec1[i] );
		len2 += double( vec2[i] * vec2[i] );
	}
	return T( inner_product / ( sqrt( len1 ) * sqrt( len2 ) ) );
}



#endif // !VECTOR_OPERATIONS_H
