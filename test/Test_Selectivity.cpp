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

using namespace cstar::selectivity::sizebased;

/*
 * Simple testing routine to confirm that functions defined within 'selectivity'.hpp have been implemented properly.
 */

int main(int argc, char** argv) {
    var output;
    var a = rand();
    var b = rand();
    var c = rand();
    var d = rand();
    var e = rand();

    output = Selectivity_1<var > (a, b, c);
    std::cout << "Variable A = " << a << "\n";
    std::cout << "Variable B = " << b << "\n";
    std::cout << "Variable C = " << c << "\n";
    std::cout << "Selectivity = " << output << "\n";
 
    return 0;
}