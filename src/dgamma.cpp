/**
*
* \file dgamma.cpp
* \brief Gamma density function
* \ingroup CSTAR
*
* \author Athol Whitten & Steve Martell
* \date 10/25/2013
*
 */

#include "../include/cstar.h"

// =========================================================================================================

dvar_vector dgamma(const dvector& x, const dvariable& a, const dvariable& b)
  {
    //returns the gamma density with a & b as parameters
    RETURN_ARRAYS_INCREMENT();
    dvariable t1 = 1./(pow(b,a)*mfexp(gammln(a)));
    dvar_vector t2 = (a-1.)*log(x)-x/b;
    RETURN_ARRAYS_DECREMENT();
    return(t1*mfexp(t2));
  }

dvariable dgamma(const prevariable& x, const double& a, const double& b)
  {
    //returns the gamma density with a & b as parameters
    RETURN_ARRAYS_INCREMENT();
    dvariable t1 = 1./(pow(b,a)*mfexp(gammln(a)));
    dvariable t2 = (a-1.)*log(x)-x/b;
    RETURN_ARRAYS_DECREMENT();
    return(t1*mfexp(t2));
  }

// =========================================================================================================
