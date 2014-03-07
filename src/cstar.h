/**
*
* \file cstar.h
* \brief Common Stock Assessment Routines
* \ingroup CSTAR
* \remarks
* This library contains declarations for a variety of
* statistical functions commonly used in stock assessment. 
* They can be used in ADMB TPL files.
* 
* \author Athol Whitten
* 
* \date 03/01/2014
*
*/

#ifndef CSTAR_H
#define CSTAR_H

#include <admodel.h>
#include "generic.cpp"
#include "dpois.cpp"
#include "dgamma.cpp"
#include "dnbinom.cpp"
#include "dmultifan.cpp"
#include "approx.cpp"

// =========================================================================================================
// Generic functions: in 'generic.cpp'
// =========================================================================================================

// Get mean length for variable objects:
double mn_length(const dvector& pobs, const dvector& mlen);
double mn_length(const dvar_vector& pobs, const dvector& mlen);

// Get standard deviation of mean length vector:
double sd_length(const dvector& pobs, const dvector& len, const dvector& mlen);

// Get normalized residulas of composition data given sample size:
dvector norm_res(const dvector& pred, const dvector& obs, double m);

// Get standard deviation of normalized residuals given observed and predicted proportions:
double sd_norm_res(const dvar_vector& pred, const dvector& obs, double m);

// Get effective sample size:
double eff_N(const dvector& pobs, const dvar_vector& phat);

// Return penalised positive values for some given vector:
dvar_vector posfun(const dvar_vector& x, const double& eps, dvariable& pen);

// Return a vector of postitions of first matches of a given value x, in a table:
ivector match(const ivector& x, const ivector& table);


// =========================================================================================================
// Density functions: in 'function_name.cpp'
// =========================================================================================================

// Poisson density function:
dvariable dpois(const dvector& k, const dvar_vector& lambda);

// Gamma density function:
dvariable dgamma(const prevariable& x, const double& a, const double& b);
dvar_vector dgamma(const dvector& x, const dvariable& a, const dvariable& b);

// Multifan-style density function:
dvariable dmultifan(const dvector& o, const dvar_vector& p, const double& s);

// Negative binomial density function:
dvariable dnbinom(const dvector& x, const dvar_vector& mu, const prevariable& k);


// =========================================================================================================
// Other functions: in 'function_name.cpp'
// =========================================================================================================

// Linearly interpolate given data points:
double approx(const double& v, const dvector& x, const dvector& y);
dvariable approx(const double& v, const dvector& x, const dvar_vector& y);


// =========================================================================================================
// Selectivity functions: in 'function_name.cpp'
// =========================================================================================================

class Selex{
private:
    dvariable m_mu;
    dvariable m_sd;
        
public:
    ~Selex() {}  // Destructor
        
    Selex( const dvariable& mu = 0., const dvariable& sd = 1.0) // default constructor
    {
        m_mu = mu;
        m_sd = sd;
    }
    
    // Logistic function
    dvector     logistic( const dvector& x,const double& mu, const double& sd );
    dvar_vector logistic( const dvector& x,const dvariable& mu, const dvariable& sd );
    
    
    // Exponential logistic
    dvector     eplogis(const dvector& x, const double& x1, const double& x2, const double& gamma);
    dvar_vector eplogis(const dvar_vector& x, const dvariable& x1, const dvariable& x2, const dvariable& gamma);
    
    
    // Linear interpolation
    dvector     linapprox(const dvector& x, const dvector& y, const dvector& xout);
    dvar_vector linapprox(const dvector& x, const dvar_vector& y, const dvector& xout);
    
    dvariable GetMu()  { return m_mu; }
    dvariable GetSd()  { return m_sd; }
};

// =========================================================================================================

#endif	/* CSTAR_H */
