# Symbol Level Simulation

## Overview
This project evaluates a single-input multi-output (SIMO) system at millimeter-wave (mmWave) and terahertz (THz) frequencies. The simulation contains realistic models of RF front-end (RFFE) components for the receiver. The RFFE models have been extracted from detailed circuit-level simulations, representing an extremely realistic characterization of the true circuit performance. For every RFFE configuration we measure the end-to-end performance and the receiver power consumption. We develop a generic model for characterizing the performance using two key parameters.

## Files
The project contains the following MATLAB scripts, classes, and data:
- ```aqnmSim.m```: This is the main simulation class. For a given RFFE configuration generates random i.i.d. transmit data and channel. Then, measures the output SNR and different receive input signal power levels. It reports the effective noise figure and power consumption of the receiver.
- ```aqnmTest.m```: This is the main simulation script. For each RFFE configuration creates a simulation object. Based on the data generated by the simulation, creates a mathematical model for characterizing the performance. The simulation data and model are saved in *.mat* file.
- ```aqnmPlot.m```: This script generates plots for analyzing the components, the end-to-end performance and the power consumption for each RFFE configuration.
- ```MultiInput.m```: This class replicates a MATLAB system object with multiple inputs.
- ```rffe_140GHz.mat```: This file contains the RFFE models for 140 GHz.

## Dependencies
This project is developed and tested on MATLAB R2020a. The project depends on the following toolbox and packages:
- Simulink
- Communications Toolbox
- Statistics and Machine Learning Toolbox
- mmwComm Package

## Acknowledgements
The detailed simulation models of the RFFE components have been developed by:
- Navid Hosseinzadeh, University of California Santa Barbara (UCSB)
- James F. Buckwalter, University of California Santa Barbara (UCSB)