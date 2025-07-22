#include <math.h>

void c_add( double *a1, double *b, double *r )
{
	*r = *a1 + *b;
	*( r + 1 ) = *( a1 + 1 ) + *( b + 1 );
}

void c_subtract( double *a1, double *b, double *r )
{
	*r = *a1 - *b;
	*( r + 1 ) = *( a1 + 1 ) - *( b + 1 );
}

void c_multiply( double *a1, double *b, double *r )
{
	*r = *a1 * *b - *( a1 + 1 ) * *( b + 1 );
	*( r + 1 ) = *a1 * *( b + 1 ) + *( a1 + 1 ) * *b;
}

void a_copy( double *s1, double *s2, int i )
{
	while( i > 0 ){
	  i--;
	  *( s2 + i ) = *( s1 + i );
	  }
}

void qcf2( double *s, int p )
{
	double s1[ 2 ], s2[ 2 ];

	c_add( s, s + 2 * p, s1 );
	c_subtract( s, s + 2 * p, s2 );
	a_copy( s1, s, 2 );
	a_copy( s2, s + 2 * p, 2 );
}

int inv_nom( char m, int x )
{
	int y = 0;

	while( m > 0 ){
	  m--;
	  y = y << 1;
	  if( x & 1 )
		y++;
	  x = x >> 1;
	  }
	return y;
}

void a_exchg( double *s0, double *s1 )
{
        double f;

        f = *s0;
        *s0 = *s1;
        *s1 = f;
        f = *(s0+1);
        *(s0+1) = *(s1+1);
        *(s1+1) = f;
}

void n_sort( double *s, int maxN )
{
	int i, k, m;

	i = maxN;
	m = 1;
	while( ( i = i >> 1 ) ^ 1 )
	  m++;
	for( i = 0; i < maxN; i++ ){
	  k = inv_nom( m, i );
          if( k > i )
	    a_exchg( s + 2 * i, s + 2 * k );
	  }
}

void fft( double *s, int maxN )
{
	int i, k, l, p;
	double r[ 2 ], v[ 2 ];
	double const Pi = 3.141593;

        n_sort( s, maxN );
	for( k = 2; k <= maxN; ){
	  p = k / 2;
	  for( i = 0; i < maxN / k; i++ )
		for( l = 0; l < p; l++ )
		  qcf2( s + 2 * ( l + i * k ), p );
	  k *= 2;
	  if( k > maxN )
		break;
	  for( i = 0; i < k / 2; i++ )
		for( l = 0; l < maxN / k; l++ ){
		  v[ 0 ] = cos( 2 * Pi * i / k );
		  v[ 1 ] =-sin( 2 * Pi * i / k );
		  c_multiply( s+2*k*l+k+2*i, v, r );
		  a_copy( r, s+2*k*l+k+2*i, 2 );
		  }
	  }
}

double length_s( double *s )
{
	return sqrt( *s * *s + *(s+1) * *(s+1) );
}
