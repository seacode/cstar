/**
*
* \file dmultifan.cpp
* \brief Multifan-style density function
* \ingroup CSTAR
*
* \author Athol Whitten & Steve Martell
* \date 10/25/2013
*
 */

#include "cstar.h"

// =========================================================================================================
  
dvariable dmultifan(const dvector& o, const dvar_vector& p, const double& s)
{
    /*
    o is the observed numbers at length
    p is the predicted numbers at length
    s is the minimum sample size
    */
    RETURN_ARRAYS_INCREMENT();      
    int lb     = o.indexmin();
    int nb     = o.indexmax();
    int I      = (nb-lb)+1;
    double n   = sum(o);
    if(min(n,s)<=0)
    {
        RETURN_ARRAYS_DECREMENT();
        return(0);
    } 
    double tau = 1./min(n,s);

    dvariable nll;
    dvector O     = o/sum(o);
    dvar_vector P = p/sum(p);
    dvar_vector epsilon(lb,nb);
    epsilon       = elem_prod(1.-P,P);

    dvariable T1,T2,T3;
    T1 = -0.5 * sum(log( 2.*M_PI*(epsilon+0.1/I) ));
    T2 = -0.5 * I * log(tau);
    T3 = sum( log(exp(-1.0 * elem_div(square(O-P),2.0*tau*(epsilon+0.1/I)) )+0.01) );
    
    nll = -1.0*(T1 + T2 + T3);
    
    //cout<<"dmultifan like = "<<nll<<"\tn = "<<1./tau<<"\tsum(P) = "<<sum(O)<<endl;
    //if(n==112) cout<<p<<endl;
    RETURN_ARRAYS_DECREMENT();
    return nll;
}

// =========================================================================================================
