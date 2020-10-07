# 3GPP NR Multi-Carrier Link-Layer Simulation

## Overview

This project evaluates a practical mobile receiver at millimeter-wave (mmWave) and terahertz (THz) frequencies. We consider a downlink system with a single NR basestation (gNB), and a single UE device. Although the NR standard was developed for 5G, the system is flexible and provides a natural baseline design for 6G systems as well. The simulation contains realistic models of RF front-end (RFFE) components for the receiver. The RFFE models have been extracted from detailed circuit-level simulations, representing an extremely realistic characterization of the true circuit performance. For each parameter setting we create a UE device that will process the received data at several input power levels.

The gNB can either use the entire wideband bandwidth by performing carrier aggregation, or transmit a single component carrier. We consider a single path channel between the gNB and the UE with random gain and phase. The gNB generates for each component carrier a physical downlink shared channel (PDSCH) that includes the information data and some physical layer signals. The demodulation reference signals (DM-RS) and the phase tracking reference signals (PT-RS). To compensate for the common phase error (CPE), 3GPP 5G NR introduced PT-RS. The UE uses DM-RS and PT-RS to perform practical channel estimation. For each slot we generate a random angle-of-departure (AoD) and angle-of-arrival (AoA) in both azimuth and elevation angles.

## Files
The project contains the following MATLAB scripts, classes, and data:
- ```NRTest.m```: This is the main simulation script.
- ```NRgNBTx.m```: This is a 3GPP NR gNB transmitter class.
- ```NRUERx.m```: This is a 3GPP NR UE receiver class.
- ```PDSCHSimParam.m```: This class contains the parameters of the PDSCH, carrier and waveform configuration of the gNB and the UE.
- ```MultiInput.m```: This class replicates a MATLAB system object with multiple inputs.
- ```rffe_140GHz.mat```: This file contains the RFFE models for 140 GHz.

## Dependencies
This project is developed and tested on MATLAB R2020a. The project depends on the following toolbox and packages:
- Simulink
- 5G Toolbox
- Antenna Toolbox
- DSP System Toolbox
- Fixed-Point Designer
- Communications Toolbox
- Phased Array System Toolbox
- Signal Processing Toolbox
- mmwComm Package

## Acknowledgements
The detailed simulation models of the RFFE components have been developed by:
- Navid Hosseinzadeh,  University of California Santa Barbara (UCSB)
