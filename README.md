# Cstar Version 1.1

by Athol Whitten, March 2014

## The Project

The Cstar (Common Stock Assessment Routines) project aims to collect, collate, and document common fisheries stock assessment modeling routines and to make them available as an open-source function library for ADMB (Auto-differentiation Model Builder). ADMB is a powerful software package for the development of nonlinear models and is freely available from the [ADMB-project website](http://www.admb-project.org). 

The library should simplify the development of fisheries stock assessment models and fish population models of all varieties. Development will be focussed on commonly used size-based functions for the assessment of hard-to-age species. For more information, please contact [Athol Whitten](mailto:whittena@uw.edu) at the School of Aquatic and Fishery Sciences, University of Washington.

The Cstar library can be used for building new stock-specific population and assessment models, or for simplifying existing ones. As an example, Cstar functions have been used, in part, to build a Generalized Model for Alaskan Crab Stocks (Gmacs). More information on Gmacs can be found at its [Github repository](https://github.com/awhitten/gmacs). The Gmacs source code includes calls to pre-existing functions from the Cstar library, eliminating the need for the main code to contain extra lines that define commonly used functions.


## Collaborators ##

Interested persons can collaborate on this open-source project by contacting [Athol Whitten](mailto:whittena@uw.edu). Collaborators may work on their own branch of the Cstar repository and make a pull request when comitting changes to the code. Some collaboators may wish to provide their Github username and be given read-write access to the repository.

### Cstar Coding Standards

1) Assume all base types are double

2) Use Google style guide: http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml

## Users ##

The entire Cstar repository can be downloaded to a local machine by clicking the `Download Zip` button on the main repository page. Github users can clone the repository using their preferred method. 

Instructions for use:

* Download the Cstar repository and save it to your local machine.

* Copy the `Cstar` folder to a directory that is in your PATH or can be referred to in your ADMB code, and use the `#include<cstar.h>` pre processor command in the GLOBALS_SECTION of your TPL file.

* Cstar functions should now work from directly within your model code.