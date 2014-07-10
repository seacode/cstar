/**
*
* \file dpois.cpp
* \brief Poisson density function
* \ingroup CSTAR
*
* \author Athol Whitten & Steve Martell
* \date 10/25/2013
*
 */

#include "../include/cstar.h"

// =========================================================================================================

dvariable dpois(const dvector& k, const dvar_vector& lambda)
{
    RETURN_ARRAYS_INCREMENT();
    int i;
    int n = size_count(k);
    dvariable nll=0;
    for(i = 1; i <= n; i++)
    {
        nll -= k(i)*log(lambda(i))+lambda(i)+gammln(k(i)+1.);
    }
    RETURN_ARRAYS_DECREMENT();
    return nll;
}

// =========================================================================================================
