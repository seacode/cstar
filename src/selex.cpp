/**
*
* \file selex.cpp
* \brief Various selectivity functions 
* \ingroup CSTAR
*
* \author Athol Whitten & Steve Martell
* \date 01/20/2014
*
 */

#include <admodel.h>

// =========================================================================================================
// Logistic function (basic with Mean and SD)                                                                  //
// =========================================================================================================

dvar_vector Selex::logistic( const dvector& x, const dvariable& mu, const dvariable& sd )
{
    return 1./(1.+mfexp(-(x-mu)/sd) );
}

dvector Selex::logistic( const dvector& x, const double& mu, const double& sd )
{
    return 1./(1.+mfexp(-(x-mu)/sd) );
}

// =========================================================================================================
// Logistic function (parameterized by size at 5 and 95% selectivity)                                                                  //
// =========================================================================================================
    
slope_par = -log(19)* selex_parms(ipnt +2);
seltmp    = 1.0 / (1.0 + mfexp(slope_par) * (length-selex_parms(ipnt +1)));
temp      = seltmp(nclass);
seltmp    /= temp;

// =========================================================================================================
// Parameter-per-size-class function                                                                    //
// =========================================================================================================

dvar_vector Selex::parmpsc(const dvar_vector& x, const dvar_vector p)
// above needs the vector of lengths  (seltmp(1,nclass), and the nclass number of parameters)

for (int iclass=1; iclass<=nclass; iclass++)
seltmp(iclass) = 1.0 / (1.0+mfexp(selex_parms(ipnt + iclass)));
temp   = seltmp(nclass);
seltmp /= temp;

// =========================================================================================================
// Exponential Logistic                                                                 //
// =========================================================================================================

dvar_vector Selex::eplogis(const dvar_vector& x, const dvariable& x1, 
                            const dvariable& x2, const dvariable& gamma)
{
    //exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
    /*
    A modified version of the exponential logistic presented in Thompson's 1994 paper in CJFAS
    
    Here the arguemnts x1 and x2 represnt the inflection points for the ascending 
    and desciending limb of the selectivity function.  Gamma is descending limb 
    parameter where gamma=0 is logistic, and gamma <1.0 dome=shaped and gamma==1 is undefined.
    
    */
    RETURN_ARRAYS_INCREMENT();
    dvariable k1,k2,alpha,beta;
    dvariable t1 = 2.-4.*gamma+2.*gamma*gamma;
    dvariable t3 = 1.+2.*gamma-2*gamma*gamma;
    dvariable t5 = sqrt(1. + 4.*gamma - 4.*gamma*gamma);
    
    k1    = log(t1/(t3+t5));
    k2    = log(t1/(t3-t5));
    beta  = (k1*x2-x1*k2)/(k1-k2);
    alpha = k2/(x2-beta);
    dvar_vector sx = (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
    
    RETURN_ARRAYS_DECREMENT();
    return sx;
}

dvector Selex::eplogis(const dvector& x, const double& x1, 
                        const double& x2, const double& gamma)
{
    //exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
    /*
    A modified version of the exponential logistic presented in Thompson's 1994 paper in CJFAS
    
    Here the arguemnts x1 and x2 represnt the inflection points for the ascending 
    and desciending limb of the selectivity function.  Gamma is descending limb 
    parameter where gamma=0 is logistic, and gamma <1.0 dome=shaped and gamma==1 is undefined.
    
    */
    double k1,k2,alpha,beta;
    double t1 = 2.-4.*gamma+2.*gamma*gamma;
    double t3 = 1.+2.*gamma-2*gamma*gamma;
    double t5 = sqrt(1. + 4.*gamma - 4.*gamma*gamma);
    
    k1    = log(t1/(t3+t5));
    k2    = log(t1/(t3-t5));
    beta  = (k1*x2-x1*k2)/(k1-k2);
    alpha = k2/(x2-beta);
    dvector sx = (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
    
    return sx;
}


// =========================================================================================================
// Linear Interpolation using approx function from R libraries ('approx1 functions defined in generic.hpp')                         //
// =========================================================================================================

dvector Selex::linapprox(const dvector& x, const dvector& y, const dvector& xout)
{
    // Piece-wise linear approximation for n points in xout between min(x) and max(y):
    int i,j,k;
    int n = xout.indexmax() - xout.indexmin() + 1;
    double v;
    dvector yout(xout.indexmin(),xout.indexmax());
    
    for(k = xout.indexmin(); k <= xout.indexmax(); k++)
    {
        v       = xout(k);
        yout(k) = approx1(v,x,y);
    }
    return yout;
}

dvar_vector Selex::linapprox(const dvector& x, const dvar_vector& y, const dvector& xout)
{
    // Piece-wise linear approximation for n points in xout between min(x) and max(y):
    int i,j,k;
    int n = xout.indexmax() - xout.indexmin() + 1;
    double v;
    dvar_vector yout(xout.indexmin(),xout.indexmax());
    for(k = xout.indexmin(); k <= xout.indexmax(); k++)
    {
        v       = xout(k);
        yout(k) = approx1(v,x,y);
    }
    return yout;
}


#endif	/* SELECTIVITY_HPP */
