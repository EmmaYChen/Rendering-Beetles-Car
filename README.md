CS348 Final Project： Glittering Beetles Car
===============
Author: Yao Chen, Xinru Hua


Introduction:
--------------
We want to render an image that contains objects like snow or beetle toy which exhibit small-scale phenomena observed as bright sparkling or glittering surface features.

We added two class, in total 5 files in pbrt, a new material class in PBRT called Flake, and a new BXDF class called FlakeBxDF. 
You could final them in src/core/material folder:
* flake.h
* flake.cpp
* glitter.h
* glitter.cpp
* sampleHSL.h

Project result:
--------------------

1. Folder 'Final Scene' : Final scene of a glittering beetles car
2. Image'final scene' : Final image for the Image Synthesis 
3. Folder 'results' : Experimental results of glittery surface
4. FInal report

How to use it
-------------
In the command line tool
* run : cmake  [path-to-pbrt]
* make: make j8




