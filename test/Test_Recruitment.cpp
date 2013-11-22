/* 
 * File:    TestSelex.cpp
 * Authors: Athol Whitten & Mathew Supernaw
 *
 * Created on July 3, 2013, 2:30 PM
 */

#include <iostream>
#include <cstdlib>

#include "Cstar.hpp"

using namespace std;

typedef double var;

using namespace cstar::recruitment::agebased;

/*
 * Simple testing routine to confirm that functions defined within 'selectivity'.hpp have been implemented properly.
 */

int main(int argc, char** argv) {
    var output;
    var a = rand();
    var b = rand();
    var c = rand();
    var d = rand();

    output = cstar::recruitment::agebased::Ricker<var > (a, b, c, d);
    std::cout << output << "\n";
    output = cstar::recruitment::agebased::BevertonHolt<var > (a, b, c, d);
    std::cout << output << "\n";
    output = cstar::recruitment::agebased::BevertonHoltConstrained<var > (a, b, c, d);
    std::cout << output << "\n";
    output = cstar::recruitment::agebased::HockeyStick<var > (a, b, c, d, d);
    std::cout << output << "\n";
    output = cstar::recruitment::agebased::Survival<var > (a);
    std::cout << output << "\n";

    return 0;
}