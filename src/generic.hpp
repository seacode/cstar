/* 
 * File:    Generic.hpp
 * Authors: Athol Whitten & Mathew Supernaw
 *
 * Created on July 3, 2013, 2:30 PM
 */

#ifndef GENERIC_HPP
#define	GENERIC_HPP

#include <cmath>

namespace cstar {

    template<class T>
    T mfexp(const T & x) {
        T b = T(60);
        if (x <= b && x >= T(-1) * b) {
            return std::exp(x);
        } else if (x > b) {
            return std::exp(b)*(T(1.) + T(2.) * (x - b)) / (T(1.) + x - b);
        } else {
            return std::exp(T(-1) * b)*(T(1.) - x - b) / (T(1.) + T(2.) * (T(-1) * x - b));
        }
    }
}


#endif	/* GENERIC_HPP */