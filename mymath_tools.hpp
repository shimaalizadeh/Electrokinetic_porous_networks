#ifndef MYMATH_TOOLS_HPP
#define MYMATH_TOOLS_HPP

#include <cmath>

double compute_norm(double *x, int n);
double heaviside(double x);
double sign(double x);

double compute_norm(double *x, int n){

	double l2norm = 0.;
	for(int i = 0; i < n; i++){
		l2norm += x[i]*x[i];
	}

	l2norm = sqrt(l2norm/n);

	return l2norm;
}

inline double heaviside(double x){
  double y;

  if(x>0)
    y=1.;

  else
    y=0.;

  return(y);
}

inline double sign(double x){
  int y;

  if(x>0 || x==0)
    y=1;

  else
    y=-1;

  return(y);
}

#endif

