/* 
 * File:    Selectivity.hpp
 * Authors: Athol Whitten & Mathew Supernaw
 *
 * Created on July 3, 2013, 2:30 PM
 */


#ifndef SELECTIVITY_HPP
#define	SELECTIVITY_HPP


/************************************
 * Common Stock Assessment Routines *
 ************************************/

namespace cstar {

    namespace selectivity {

        namespace sizebased {

                template<class T,class TT>
                T Selectivity_1(const T &a, TT length, const T &selex_parameter) {
                   
                    return 1.0 / (1.0 + mfexp(a * (length - selex_parameter)));
                }
                
                 template<class T>
                T Selectivity_2(const T &a, const T &selex_parameter) {
                   
                    return 1.0 / (1.0 + mfexp(a *  selex_parameter));
                }
            
            } //end sizebased
        
        } // selectivity
    
    } // end cstar

#endif	/* SELECTIVITY_HPP */

