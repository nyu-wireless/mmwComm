classdef NRUERx < matlab.System
    % 5G NR UE receiver class
	
	properties                
        carrierConfig;	% carrier configuration         
        pdschConfig;	% default PDSCH configuration
        waveformConfig;	% waveform configuration
        ofdmGrid;		% before pre-coding nsc x nsym x nlayers
        pdschSymTx;     % Tx symbols
		
		% Received symbols
        pdschChanEst;   % Channel estimate on the PDSCH
        pdschSymRaw;    % Raw symbols before equalization
        pdschSymEq;     % Equalized symbols
		
		% Channel and noise estimate
        noiseEst;
		wrx;
        Nslot = 0;	% slot number
		nrx;		% number of rx antennas
		
		% RFFE
		rffe;		% rf front-end object
		isLinear;	% include the rffe non-linearities
		lnaGain;	% lna gain in dB
		lnaNF;		% lna noise figure in dB
		lnaAmpLut;	% lna fund tone power curve
		lnaPower;	% lna power in mW
		
		mixGain;	% mixer gain in dB
		mixNF;		% mixer noise figure in dB
		mixAmpLut;	% mixer fund tone power curve
		mixPLO;		% mixer input power from the local oscillator
		mixPower;	% mixer power in mW
				
		% ADC
		fsamp;				% sample frequency of the ADC
		adcFOM = 65e-12;	% FOM of the ADC 65 fJ/conv
		nbits = 4;			% ADC bits per dimension
		
		% PHY
		phy;			% Physical layer object
		isFixPoint;		% use fixed-point arithmetic
		enableCA;		% enable carrier aggregation
        ncc;			% number of component carriers
		
		snrInTest = 0;	% input SNR test vector
		debug = false;	% 'true' to display debug messages
	end
	
	methods
		function obj = NRUERx(varargin)
			% Constructor
			
			% Set parameters from constructor arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
			
			% Create the RX RFFE
			obj.rffe = mmwsim.rffe.RFFERx(...
				'nrx', obj.nrx, ...
				'lnaNF', obj.lnaNF, ...
				'lnaGain', obj.lnaGain, ...
				'lnaPower', obj.lnaPower, ...
				'lnaAmpLut', obj.lnaAmpLut, ...
				'mixNF', obj.mixNF, ...
				'mixPLO', obj.mixPLO, ...
				'mixGain', obj.mixGain, ...
				'mixPower', obj.mixPower, ...
				'mixAmpLut', obj.mixAmpLut, ...
				'fsamp', obj.fsamp, ...
				'nbits', obj.nbits, ...
				'isLinear', obj.isLinear);
		end
		
		function NF = nf(obj)
			% Calculate the effective noise figure
			NF = obj.rffe.nf();
		end
		
		function P = power(obj)
			P = obj.rffe.power();
		end
	end
	
	methods (Access = protected)
		
		function setupImpl(obj)
			
 			% Create the RX PHY			
			obj.phy = mmwsim.phy.PHYRx(...
  				'ncc', obj.ncc, ...
  				'isFixPoint', obj.isFixPoint, ...
				'enableCA', obj.enableCA, ...
  				'NRB', obj.carrierConfig.NRB, ...
  				'SCS', obj.carrierConfig.SubcarrierSpacing, ...
  				'fsamp', obj.waveformConfig.SamplingRate, ....
				'nbadc', obj.nbits);
		end
		
		function [snrOut] = stepImpl(obj, y)
            % Find the spatial signatures and element gains of each path 
			% based on their angles of arrival.
 			nsnr = length(obj.snrInTest);
 			snrOut = zeros(nsnr, obj.ncc);
 			
 			for isnr = 1:nsnr
 				% Get the SNR and scale the input signal
 				ynoisy = 10^(0.05*obj.snrInTest(isnr))*y;
 		
 				% Pass the data from the RFFE
 				ynoisy = obj.rffe.step(ynoisy);
 
 				% Pass the data from the PHY
				ynoisy = obj.phy.step(ynoisy);

				% Get information for PDSCH, DM-RS and PT-RS allocations
				[pdschIndices,dmrsIndices,dmrsSym,ptrsIndices, ...
					ptrsSym, pdschIndicesInfo] = ...
					mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
				
				% Loop over Component Carrier
				for icc = 1:obj.ncc

					% Perform RX beamforming by multiplying y with the
					% the RX BF vector.
					z = ynoisy(:,:,icc)*obj.wrx;

					% Demodulate the RX signal
					obj.ofdmGrid = mmwsim.nr.hOFDMDemodulate(...
						obj.carrierConfig, z);

					% Get channel estimate.
					% This is a poor channel estimate since we have not done
					% carrier and timing estimation.  But, this is OK for now.
					[chanEstGrid, obj.noiseEst] = nrChannelEstimate(...
						obj.ofdmGrid, dmrsIndices, dmrsSym,...
						'CyclicPrefix', obj.carrierConfig.CyclicPrefix,...
						'CDMLengths', pdschIndicesInfo.CDMLengths);
					
					% Initialize temporary grid to store equalized symbols
					tempGrid = zeros(...
						obj.waveformConfig.NSubcarriers,...
						obj.waveformConfig.SymbolsPerSlot, ...
						obj.pdschConfig.NLayers);
					
					% Equalize the PDSCH symbols and map them to tempGrid
					[pdschRx, pdschHest] = nrExtractResources(...
						pdschIndices, obj.ofdmGrid, chanEstGrid);
					pdschEq = nrEqualizeMMSE(pdschRx, pdschHest, obj.noiseEst);
					tempGrid(pdschIndices) = pdschEq;	
					
					% Equalize the PTRS symbols and map them to tempGrid
					[ptrsRx, ptrsHest] = nrExtractResources(...
						ptrsIndices, obj.ofdmGrid, chanEstGrid);
					ptrsEq = nrEqualizeMMSE(ptrsRx, ptrsHest, obj.noiseEst);
					tempGrid(ptrsIndices) = ptrsEq;
					
					% Estimate the residual channel at the PTRS locations
					cpe = nrChannelEstimate(tempGrid, ptrsIndices, ptrsSym);
					
					% Sum estimates across subcarriers, receive antennas, and
					% layers. Then, get the CPE by taking the angle of the
					% resultant sum.
					cpe = angle(sum(cpe,[1 3 4]));
					
					% Correct CPE in each OFDM symbol within the range of 
					% reference PTRS OFDM symbols
					symLoc = pdschIndicesInfo.PTRSSymbolSet(1)+1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
					tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1j*cpe(symLoc));
					
					% Extract raw symbols and channel estimate on PDSCH
					obj.pdschSymRaw(:,icc) = obj.ofdmGrid(pdschIndices);
					obj.pdschChanEst(:,icc) = chanEstGrid(pdschIndices);
					obj.pdschSymEq(:,icc) = tempGrid(pdschIndices);
					
					% Use a linear receiver to estimate the output SNR
					%
					% xhat = a*x + d,  d ~ CN(0, E|xhat-x|^2)
					xhat = obj.pdschSymEq(:,icc);
					x = obj.pdschSymTx(:,icc);
					xvar = mean(abs(x).^2, 'all');
					a = mean(conj(xhat).*x)/mean(abs(x).^2);
					dvar = mean(abs(xhat - a*x).^2);

					% Measure the output SNR
					snrOut(isnr,icc) = 10*log10(abs(a).^2*xvar/dvar);
				end
				
				% Display progress. 
				if obj.debug
					fprintf(1, 'snrIn = %12.4e, snrOut = %12.4e\n', ...
						obj.snrInTest(isnr), mean(snrOut(isnr),2));
				end
			end
			snrOut = mean(snrOut, 2);
		end
	end
end

