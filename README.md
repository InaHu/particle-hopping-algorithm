The Particle-Hopping-Algorithm
============

This 'Particle-Hopping-Algorithm' solves the PDE-system presented in the paper "On the Role of Vesicle Transport in Neurite Growth: Modelling and Experiments". A preprint is available at https://arxiv.org/pdf/1908.02055.pdf
It plots both the time evolution of anterograde and retrograde moving particles and the time evolution of the concentration in the pools.
For an explaination see chapter 5 "Numerical Simulations" of the paper which explains the numerics in more detail.

It is advisible (but not nescessary) to use the helperFiles package provided by Janic Föcke (FAU Erlangen-Nürnberg) for a better visualization of the images provided by this code. 
Without this toolbox the images might look very ugly when you try to print them. If you do not want to use this toolbox simply uncomment the lines that use the 'exportFigure'-commands.
The helper-Files packages is free accessible under https://repo.mi.uni-erlangen.de/imaging/code/public/helperFiles (see its Read.me for explanation how to use it).


Dependencies:
-------------
* Matlab >= R2017b
* helper files package (not nescessary, see above description)


Usage:
------
This algorithm does not require a specific spot to be stored. Please choose the folder at your pleasure. To use it simply use `addpath` with the chosen save spot. 

Then all functions are available and you can execute the main-file.


Authors:
--------
* Ina Humpert([ina.humpert@uni-muenster.de](mailto:ina.humpert@uni-muenster.de))
* Jan-Frederik Pietschmann([jfpietschmann@math.tu-chemnitz.de](mailto:jfpietschmann@math.tu-chemnitz.de))

Acknowledgements:
--------
* Andreas Püschel
* Danila Di Meo
* Janic Föcke 
