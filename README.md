# PCICEFEM-Von-Karman
*Von karman vortex street by flow passing over a cylinder (UNDER DEVELOPMENT)*

This repository contains code that solves for the 2D flow of air around a cylinder geometry in an attempt to create a wake of air called the Von Karman Vortex street. The code is as yet incomplete , and riddled with bugs, but is executable. Explained here in this readme are some important files.

# von karmann.igs and von karmann.neu
The .igs and .neu files contain the grid of the domain.

# flow.in
Contains inflow conditions: Pressure, Temperature, Ratio of specific heats, molecular weight, prandtl number, mach number, viscosity, CFL number and restart

# inflow.dat, outflow.dat, wall.dat
Contains boundary nodes with the respective properties.

# restart.in
Save file for all parameters, in case simulations need to be carried on further from the point where it was saved. Simulation is restarted by setting the restart number to 1 in flow.in

# FENSES.CPP
This is the master file that houses the solver code. All input files, and boundary data files are read, the neighboring nodes to each node is identified, properties of elements such as area and shortest edge are deduced. Then the problem is initialized by assessing the initial values of all parameters such as x and y velocities, pressure, temperature, energy and density. The boundary conditions are put in place next, having identified how each kind of boundary nodes are going to be solved for, then the timestep is calculated and the first set of equations formed (in matrices), before finally starting the time-step iterations. Per time-step iteration, the values are solved for implicitly iteratively, to increase accuracy and maintain overall computational stability. Finally the results are written onto results.dat, restart.dat (in case this is from where iterations are to be run again) and tec.dat.

# tec.dat
This file is read in tecplot (or paraview) to visualize the results. Temperature contours.png and Mach contours.png are the respective visualizations of the domain.
