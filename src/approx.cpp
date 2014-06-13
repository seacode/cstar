/**
*
* \file approx.cpp
* \brief Linearly interpolate given data points
* \ingroup CSTAR
*
*  Return a list of points which linearly interpolate 
*  given data points, or a function performing the 
*  linear (or constant) interpolation [from R functions].
*
* \author Steve Martell & Athol Whitten
* \date 12/15/2013
*
 */

#ifndef APPROX_H
#define APPROX_H
#include "../include/cstar.h"
// =========================================================================================================
  
double approx(const double& v, const dvector& x, const dvector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    i = x.indexmin();
    j = x.indexmax() - 1;

    /* handle out-of-domain points */
    if(v < x[i]) return min(y);
    if(v > x[j]) return max(y);

    /* find the correct interval by bisection */
    while(i < j - 1) 
    { /* x[i] <= v <= x[j] */
        ij = (i + j)/2; /* i+1 <= ij <= j-1 */
        if(v < x[ij]) j = ij;
        else i = ij;
        /* still i < j */
    }
    /* provably have i == j-1 */

    /* interpolation */
    if(v == x[j]) return y[j];
    if(v == x[i]) return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

    /* linear */
    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
}/* approx() */


dvariable approx(const double& v, const dvector& x, const dvar_vector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    i = x.indexmin();
    j = x.indexmax() - 1;

    /* handle out-of-domain points */
    if(v < x[i]) return min(y);
    if(v > x[j]) return max(y);

    /* find the correct interval by bisection */
    while(i < j - 1) 
    { /* x[i] <= v <= x[j] */
        ij = (i + j)/2; /* i+1 <= ij <= j-1 */
        if(v < x[ij]) j = ij;
        else i = ij;
        /* still i < j */
    }
    /* probably have i == j-1 */

    /* interpolation */
    if(v == x[j]) return y[j];
    if(v == x[i]) return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

    /* linear */
    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
}/* approx() */
#endif
// =========================================================================================================
  