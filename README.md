# Millimeter Wave Communications Package

## Overview

The goal of `mmwComm` is to provide a single platform for detailed link-layer simulation and FPGA design based in MATLAB and [Xilinx System Generator](https://www.mathworks.com/products/connections/product_detail/xilinx-system-generator-for-dsp.html).  The package is intended to assist the design and evaluation of various mmWave systems including:
   * Millimeter wave communication and sensing systems 
   * Systems in 5G bands (e.g. 28, 37 and 73 GHz) as well as systems above 100 GHz
   * MIMO communication front-end systems
   * Standards compliant systems including 5G NR or 802.11-type systems
   * MIMO channel sounding 
   * MIMO radar   
   
The eventual desired features are:
* Detailed simulation of the RFFE including mixer, LNA and PA non-linearities and phase noise
* Support arbitrary antenna array layouts and element patterns (using the [phase array toolbox](https://www.mathworks.com/products/phased-array.html) )
* Simulation and design of the key PHY components including digital filtering, OFDM processing, equalization, MIMO processing, channel encoding / decoding
* Floating and bit fixed point simulation
* Key digital blocks are written in  [Xilinx System Generator](https://www.mathworks.com/products/connections/product_detail/xilinx-system-generator-for-dsp.html) to enable direct bit-exact compiling to FPGA blocks
* Support small numbers of nodes link-level simulation of multi-user and interference studies

Importantly, the package is *not* intended for:
*  MAC and upper layer simulation.  The simulation will be too slow. 
*  DSP code is not automatically generated.  Instead, any DSP code will be modeled in MATLAB, and will require manual transfer to the DSP or ARM core.  Only the FPGA code is automatically generated.