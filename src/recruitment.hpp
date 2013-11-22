/* 
 * File:    Recruitment.hpp
 * Authors: Athol Whitten & Mathew Supernaw
 *
 * Created on July 3, 2013, 2:30 PM
 */

#ifndef CSTAR_RECRUITMENT_HPP
#define	CSTAR_RECRUITMENT_HPP

#include <stdlib.h>


/************************************
 * Common Stock Assessment Routines *
 ************************************/

namespace cstar {

    namespace recruitment {

        /**
         *Recruitment functions based on age.
         */
        namespace agebased {

            /**
             * 
             * @param MinPoss
             * @param MaxPoss
             * @param Inflec
             * @param Xvar
             * @param Y1
             * @param Y2
             * @return 
             */
            template<class T>
            static T Join_Fxn(const T& MinPoss, const T& MaxPoss, const T& Inflec, const T& Xvar, const T& Y1, const T& Y2) {
                T Yresult;
                T join;
                join = CSTAR_REAL(1.000) / (CSTAR_REAL(1.000) + mfexp(CSTAR_REAL(1000.0) * (Xvar - Inflec) / (MaxPoss - MinPoss))); //  steep joiner at the inflection
                Yresult = Y1 * (join) + Y2 * (CSTAR_REAL(1.000) - join);
                return Yresult;
            }

            /**
             * Implementation of the Previous Placement recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T PreviousPlacement(const T& SPB_current) {
                std::cout << "Critical error:  see warning" << std::endl;
//#warning will eventually add standard logging...
                std::cerr << "B-H constrained curve is now Spawn-Recr option #6" << std::endl;
                exit(1);
            }

            /**
             * Implementation of the Ricker recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T Ricker(const T& SPB_current, const T &Recr_virgin, const T &SPB_virgin, const T &steepness) {

                T Recr_virgin_adj;
                T SPB_virgin_adj;
                T SPB_curr_adj;

                SPB_curr_adj = SPB_current + CSTAR_REAL(0.100); // robust
                Recr_virgin_adj = Recr_virgin;
                SPB_virgin_adj = SPB_virgin;
                return Recr_virgin_adj * SPB_curr_adj / SPB_virgin_adj * mfexp(steepness * (CSTAR_REAL(1.000) - SPB_curr_adj / SPB_virgin_adj));

            }

            /**
             *  Implementation of the BevertonHolt recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T BevertonHolt(const T& SPB_current, const T &Recr_virgin, const T &SPB_virgin, const T &steepness) {

                T SPB_curr_adj;
                T Recr_virgin_adj;
                T SPB_virgin_adj;
                SPB_curr_adj = SPB_current + CSTAR_REAL(0.100); // robust
                Recr_virgin_adj = Recr_virgin;
                SPB_virgin_adj = SPB_virgin;
                return (CSTAR_REAL(4.) * steepness * Recr_virgin_adj * SPB_curr_adj) / (SPB_virgin_adj * (CSTAR_REAL(1.000) - steepness)+(CSTAR_REAL(5.) * steepness - CSTAR_REAL(1.000)) * SPB_curr_adj);
            }

            /**
             *  Implementation of the None recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T None(const T& SPB_current) {

                return SPB_current + CSTAR_REAL(0.100);
            }

            /**
             * Implementation of the Hockey Stick recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T HockeyStick(const T& SPB_current, const T &Recr_virgin, const T &SPB_virgin, const T &SR_parm, const T &steepness) {

                T SPB_curr_adj;
                T Recr_virgin_adj;
                T SPB_virgin_adj;
                SPB_curr_adj = SPB_current + CSTAR_REAL(0.100); // robust
                Recr_virgin_adj = Recr_virgin;
                SPB_virgin_adj = SPB_virgin;
                //        if (SR_parm.indexmax() < 3) {
                //            cout << "Critical error:  see warning" << endl;
                //    #warning will eventually add standard logging...
                //            std::cerr << "SR_param index out of bounds." << endl;
                //            exit(1);
                //        }
                T temp = SR_parm * Recr_virgin_adj + SPB_curr_adj / (steepness * SPB_virgin_adj)*(Recr_virgin_adj - SR_parm * Recr_virgin_adj); //  linear decrease below steepness*SPB_virgin_adj
                return Join_Fxn<T > (CSTAR_REAL(0.0) * SPB_virgin_adj, SPB_virgin_adj, steepness*SPB_virgin_adj, SPB_curr_adj, temp, Recr_virgin_adj);
            }

            /**
             * Implementation of the Constrained  Beverton-Holt constrained 
             * recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T BevertonHoltConstrained(const T& SPB_current, const T &Recr_virgin, const T &SPB_virgin, const T &steepness) {


                T SPB_BH1;
                T Recr_virgin_adj;
                T SPB_virgin_adj;
                //        T steepness;
                T SPB_curr_adj;

                SPB_curr_adj = SPB_current + CSTAR_REAL(0.100); // robust
                Recr_virgin_adj = Recr_virgin;
                SPB_virgin_adj = SPB_virgin;
                ;


                if (SPB_curr_adj > SPB_virgin_adj) {
                    SPB_BH1 = SPB_virgin_adj;
                } else {
                    SPB_BH1 = SPB_curr_adj;
                }
                return (CSTAR_REAL(4.) * steepness * Recr_virgin_adj * SPB_BH1) / (SPB_virgin_adj * (CSTAR_REAL(1.) - steepness)+(CSTAR_REAL(5.) * steepness - CSTAR_REAL(1.)) * SPB_BH1);
            }

            /**
             * Implementation of the Survival recruitment model.
             * @param SPB_current
             * @return 
             */
            template<class T>
            static T Survival(const T& SPB_current) {
                std::cerr << __func__ << " Not Yet implemented." << std::endl;
                T ret = CSTAR_REAL(0);
                return ret;
            }


        }//end agebased


        /**
         *Recruitment functions based on size
         */
        namespace sizebased {


        }//end sizebased




    }//end recruitment


}//end cstar

#endif	/* RECRUITMENT_HPP */

