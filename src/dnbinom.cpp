/**
*
* \file dnbinom.cpp
* \brief Negative binomial density function
* \ingroup CSTAR
*
* \author Athol Whitten & Steve Martell
* \date 10/25/2013
*
 */

#include "../include/cstar.h"

// =========================================================================================================
  
dvariable dnbinom(const dvector& x, const dvar_vector& mu, const prevariable& k)
{
    //the observed counts are in x
    //mu is the predicted mean
    //k is the overdispersion parameter
    if (value(k)<0.0)
    {
        cerr<<"k is <=0.0 in dnbinom()";
        return(0.0);
    }
    RETURN_ARRAYS_INCREMENT();
    int i,imin,imax;
    imin=x.indexmin();
    imax=x.indexmax();
    dvariable loglike = 0.;

    for(i = imin; i<=imax; i++)
    {
        cout<<"mu "<<mu(i)<<endl;
        loglike += gammln(k+x(i))-gammln(k)-gammln(x(i)+1)+k*log(k)-k*log(mu(i)+k)+x(i)*log(mu(i))-x(i)*log(mu(i)+k);
    }
    RETURN_ARRAYS_DECREMENT();
    return(-loglike);
}

// =========================================================================================================
  