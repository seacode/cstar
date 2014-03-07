/**
*
* \file generic.cpp
* \brief Generic functions useful for model development
* \ingroup CSTAR
*
* \author Athol Whitten & Jim Ianelli
* \date 10/25/2013
*
 */

#include "cstar.h"

// =========================================================================================================

// mn_length(): Return mean length for variable objects (model estimates).

double mn_length(const dvector& pobs, const dvector& mlen)
{
  double mobs = (pobs*mlen);
  return mobs;
}

double mn_length(const dvar_vector& pobs, const dvector& mlen)
{
  double mobs = value(pobs*mlen);
  return mobs;
}

// ------------------------------------------------------------------------------------ //
// sd_length(): Return standard deviation of length.
  
double sd_length(const dvector& pobs, const dvector& len, const dvector& mlen)
{
  double mobs = (pobs*len);
  double stmp = sqrt((elem_prod(mlen,mlen)*pobs) - mobs*mobs);
  return stmp;
}

// ------------------------------------------------------------------------------------ //
// norm_res(): Returns normalized residuals of composition data given sample size.
  
dvector norm_res(const dvector& pred, const dvector& obs, double m)
{
  RETURN_ARRAYS_INCREMENT();
  pred += 0.0001;
  obs  += 0.0001;
  dvector nr(1,size_count(obs));
  nr = elem_div(obs-pred,sqrt(elem_prod(pred,(1.-pred))/m));
  RETURN_ARRAYS_DECREMENT();
  return nr;
}

// ------------------------------------------------------------------------------------ //
// sd_norm_res(): Computes standard deviation of normalized residuals given observed and predicted proportions.

double sd_norm_res(const dvar_vector& pred, const dvector& obs, double m)
{
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred)+ 0.0001;
  sdnr = std_dev(norm_res(pp,obs,m));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;
}

// ------------------------------------------------------------------------------------ //
// eff_N(): Computes effective sample size.

double eff_N(const dvector& pobs, const dvar_vector& phat)
{
  pobs += 0.0001;
  phat += 0.0001;
  dvar_vector rtmp = elem_div((pobs-phat),sqrt(elem_prod(phat,(1-phat))));
  double vtmp;
  vtmp = value(norm2(rtmp)/size_count(rtmp));
  return 1./vtmp;
}

// ------------------------------------------------------------------------------------ //
// posfun(): Return penalised positive values for some given vector.

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
// match(): Return a vector of postitions of first matches of given value x, in a table.

ivector match(const ivector& x, const ivector& table)
  {
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

// =========================================================================================================
