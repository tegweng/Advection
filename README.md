# Advection

Prerequisites:

Recommended:

Spyder (through python (x,y)

Or other python editor for 2.7.10 (Code will work for later versions of python but use notes below).

Files included in the repository:

Main file:

advection.py

Files with advection schemes:

Part I:

AdvectionFTCS.py
Advection_FTBS.py
advectionBTBS.py
advectionCTCS.py

Part II:

Artificial_diffusion.py
SemiLagrangian.py
Warming and Beam.py
advectionTVD.py

Initial conditions:

initialConditions.py

Code testing:

CodeTests.py

There is also a folder containing pdf files of the graphs for the report, a gitignore file,
this readMe file and a diagnostics.py file (not changed).

Proceedure:

Note that if you are in python 3 you should use the runfile options rather than execfile.

In order to run the advcetion codes follow the following proceedure:

1. Run the file advection.py
2. In the console, run the function main
  a. To use the initial conditions for the square wave, use main with the conditions with
  the standard conditions 
  b. To use the initial conditions for the cosine wave, use main with "func = cosine"
3. The code should produce a pdf file in the folder plots which shows the initial conditions against the advection schemes for part I and print the courant number, dx, dt, nt, end time, d_2 and d_4. 

In order to run the testing for all functions use the following proceedure:

1. Run the file 
2. All code should pass

