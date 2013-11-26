/* 
 * File:    Generic.hpp
 * Authors: Athol Whitten
 *
 * Created July 3, 2013
 * Updated November 22, 2013
 */

#ifndef GENERIC_HPP
#define	GENERIC_HPP

/************************************
 * Common Stock Assessment Routines *
 ************************************/

#include <admodel.h>

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

// ------------------------------------------------------------------------------------ //
    
dvariable dmultifan(const dvector& o,const dvar_vector& p,const double& s)
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

// ------------------------------------------------------------------------------------ //
    
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

// ------------------------------------------------------------------------------------ //

dvar_vector posfun(const dvar_vector& x, const double& eps, dvariable& pen)
  {
    int i;
    dvar_vector xp(x.indexmin(),x.indexmax());
    for(i=x.indexmin();i<=x.indexmax();i++)
    {
        if(x(i)>=eps)
        {
            xp(i) = x(i);
        }
        else
        {
            pen += 0.01*square(x(i)-eps);
            xp(i) = eps/(2.-x(i)/eps);
        }
    }
    return(xp);
  }

// ------------------------------------------------------------------------------------ //

ivector match(const ivector& x, const ivector& table)
  {
    //returns a vector of positions of first matches of x in table.
    int i,j;
    ivector pos(x.indexmin(),x.indexmax());
    for(i=x.indexmin();i<=x.indexmax();i++)
    {
        for(j=table.indexmin();j<=table.indexmax();j++)
        {
            if(x(i) == table(j) )
            {
                pos(i) = j;
                break;
            }
        }
    }
    return(pos);
  }

// ------------------------------------------------------------------------------------ //

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

// ------------------------------------------------------------------------------------ //
// Previously in 'generic.hpp', these approx functions by SM as per R libraries
// ------------------------------------------------------------------------------------ //

double approx1(const double& v, const dvector& x, const dvector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    //if(!n) return R_NaN;

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
}/* approx1() */

dvariable approx1(const double& v, const dvector& x, const dvar_vector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    //if(!n) return R_NaN;

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
}/* approx1() */

// ------------------------------------------------------------------------------------ //

#endif	/* GENERIC_HPP */