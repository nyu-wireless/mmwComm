# 3GPP NR Multi-Carrier Link-Layer Simulation

## Overview
This project evaluates a practical mobile receiver at millimeter-wave (mmWave) and terahertz (THz) frequencies. We consider a downlink system with a single NR basestation (gNB), and a single UE device. Although the NR standard was developed for 5G, the system is flexible and provides a natural baseline design for 6G systems as well. The simulation contains realistic models of RF front-end (RFFE) components for the receiver. The RFFE models have been extracted from detailed circuit-level simulations, representing an extremely realistic characterization of the true circuit performance.

We consider either a 3GPP NR multi-path channel or a single-path channel with random gain and phase between the gNB and the UE. The gNB can either use the entire wideband bandwidth by performing carrier aggregation, or transmit a single component carrier. For each component carrier generates a physical downlink shared channel (PDSCH) that includes the information data and some physical layer signals. More specifically, the demodulation reference signals (DM-RS) and the phase tracking reference signals (PT-RS). The UE uses the DM-RS and PT-RS to perform practical channel estimation. For each slot we generate a random angle-of-departure (AoD) and angle-of-arrival (AoA) in both azimuth and elevation angles. For each parameter setting we create a UE device that will process the received data at several input power levels.


## Files
The project contains the following MATLAB scripts, classes, and data:
- ```rffeTest.m```: This is the main simulation script. Using the MATLAB antenna element models creates an antenna array for the UE and the gNB. For each slot creates a gNB transmitter, a single/multi path channel, and one UE receiver for each RFFE parameter setting. The UE will process the received data at several input power levels. It creates a model for the performance using two keys parameter. The data will be saved in a *.mat* file.
- ```NRgNBTx.m```: This is a 3GPP NR gNB transmitter class. It performs carrier aggregation with ideal baseband processing. For each component carrier generates a PDSCH with DM-RS and PT-RS.
- ```NRUERx.m```: This is a 3GPP NR UE receiver class. It can have linear RFFE or non-linear RFFE with intermodulation distortion (IMD). It performs the carrier disaggregation with either ideal or fixed-point processing. It used the PT-RS to remove the common phase error (CPE) noise.
- ```PDSCHSimParam.m```: This class contains the configuration parameters of the PDSCH, carrier and waveform of the gNB and the UE.
- ```rffePlot.m```: This script generates plots for analyzing the components, the end-to-end performance and the power consumption for each RFFE configuration.
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
