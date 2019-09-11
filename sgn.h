#ifndef __SGN_H__
#define __SGN_H__
#include <assert.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <iostream>


inline double sgn(double x) 
{
    return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}

inline double mag(double x)
{
	return (x > 0.0) ? x : ((x < 0.0) ? -x : 0.0);
}

inline void Cross(double* A,double* B, double* C)
{	
	

	A[0] = B[1] * C[2] - C[1] * B[2];

    A[1] = C[0] * B[2] - B[0] * C[2];

    A[2] = B[0] * C[1] - C[0] * B[1];

	return;
}

inline void Add( double* A, double* B)
{
	
	
	A[0]+=B[0];
	A[1]+=B[1];
	A[2]+=B[2];

	return ;

}

inline void MultiplyS(double* A, double gammap)
{
	
	
	A[0]*=gammap;
	A[1]*=gammap;
	A[2]*=gammap;

	return ;

}

inline void MultiplyV( double* A, double* B)
{
	
	A[0]*=B[0];
	A[1]*=B[1];
	A[2]*=B[2];

	return ;

}

inline float FISqrt(float x) {
  float xhalf = 0.5f * x;
  int i = *(int*)&x;         // evil floating point bit level hacking
  i = 0x5f3759df - (i >> 1);  // what the fuck?
  x = *(float*)&i;
  x = x*(1.5f-(xhalf*x*x));
  return x;
}

inline void DiffMagV(double*A, double* B,double* C)
{
	

	A[0]=(B[0]-C[0])*(B[0]-C[0])+(B[1]-C[1])*(B[1]-C[1])+(B[2]-C[2])*(B[2]-C[2]);
	A[0]=A[0]*FISqrt(A[0]);
	return;

}

inline void Dott(double* out, double* A,double* B)
{
	out[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];

	return ;
}

inline void Norm(double* A)
{
	MultiplyS(A,1.0/std::sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]));
	return ;

}

inline void Normalize(double* A, double B)
{
	MultiplyS(A,B/std::sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]));
	return ;

}

#endif

