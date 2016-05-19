#ifndef MYMATH_TOOLS_H
#define MYMATH_TOOLS_H

#include <math.h>
double heaviside(double x);
double sign(double x);

inline double heaviside(double x)
{
  double y;

  if(x>0)
    y=1.;

  else
    y=0.;

  return(y);
}

inline double sign(double x)
{
  int y;

  if(x>0 || x==0)
    y=1;

  else
    y=-1;

  return(y);
}

#endif

