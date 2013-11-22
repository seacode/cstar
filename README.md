# Cstar Version 1.0 

by Athol Whitten, November 2013

## The Project

The Cstar (Common Stock Assessment Routines) project aims to collect, collate, and document common fisheries stock assessment modelling routines and to make them available as an open-source function library for C++ compatible applications, including ADMB (Auto-differentiation Model Builder). ADMB is a powerful software package for the development of nonlinear models and is freely available from the [ADMB-project website](http://www.admb-project.org). 

The library will simplify the development of fisheries stock assessment models and fish population models of all varieties. At first, development will be focussed on coding commonly used size-based functions for the assessment of hard-to-age species. For more information, please contact [Athol Whitten](mailto:whittena@uw.edu) at the School of Aquatic and Fishery Sciences, University of Washington.


## The Plan

The Cstar library can be used for building new stock-specific population and assessment models, or for simplifying existing ones. To demonstrate this utility, an example stock assessment model previously developed to assess the Bristol Bay Red King Crab Fishery (BBRKC) is provided here in its original, and simplified forms (see `Examples` folder). The simplified form is a work in progress, but includes calls to pre-existing functions from the Cstar library, eliminating the need for the main assessment code to contain extra lines that define commonly used functions.

Cstar functions can also be used to build generalized modeling packages, and will be used to build a Generalized Model for Alaskan Crab Stocks (Gmacs). More information on Gmacs can be found at its [Github repository](https://github.com/awhitten/gmacs).

## Collaborators ##

Interested persons can collaborate on this open-source project by contacting [Athol Whitten](mailto:whittena@uw.edu). Collaborators may work on their own branch of the Cstar repository and make a pull request when comitting changes to the code. Some collaboators may wish to provide their Github username and be given read-write access to the repository.

### Cstar Coding Standards

1) Assume all base types are double

2) Use template functions
        - make portable code
        - never hard-code types from some library in you functions(i.e., variable, adnumber, etc).

3) Use Google style guide: http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml

4) Always define namespaces inside cstar namespace.
  namespace cstar {
                	
                	namespace selectivity{
                	
                  }
  }

## Users ##

The entire Cstar repository can be downloaded to a local machine by clicking the `Download Zip` button on the main repository page. Github users can clone the repository using their preferred method. 

Instructions for use:

	Download the Cstar repository and save it to your local machine in a directory of your choosing. 

	Copy the `Cstar` folder to a directory that is easily referred to in your C++ code. ADMB users may wish to copy the folder to the `Contrib` folder in their ADMB directory.