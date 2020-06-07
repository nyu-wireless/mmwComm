classdef NRUERx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;  % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grid b
        ofdmGrid;     % Before pre-coding nsc x nsym x nlayers
        
        % Channel and noise estimate
        noiseEst;
        
        % Recived symbols
        pdschChanEst;   % Channel estimate on the PDSCH
        pdschSymRaw;    % Raw symbols before equalization
        pdschSymEq;     % Equalized symbols
        dmrsSym;        % DM-RS Reference symbols
        
        % Slot number
        Nslot = 0;
        
        % RX beamforming vector.  
        rxBF;
        
        % Timing offset
        offset;
	
		% RF Front-End
		rfferx;

		% Physical layer processing
		phyrx;

		% Component Carrier Aggregation
		CarrierAggregationEnable = true;
		componentCarrier = 4;
		fc = [150e6, 50e6, -50e6, -150e6];
		fsamp = 491.52e6;	% sample frequency
    end
    
    methods
        function obj = NRUERx(simParam, varargin)
            % Constructor
           
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the 
            % objects need to be converted to older structures.
            obj.carrierConfig = mmwsim.nr.objtostruct( simParam.carrierConfig );
            obj.pdschConfig = mmwsim.nr.objtostruct( simParam.pdschConfig );
            obj.waveformConfig = mmwsim.nr.objtostruct( simParam.waveformConfig );
			obj.rfferx = mmwsim.rffe.RFFERx('nbadc', simParam.nbadc);
			obj.phyrx = mmwsim.phy.PhyRx('nbadc', simParam.nbadc, ...
				'carrierConfig', simParam.carrierConfig);
            
            % Set parameters from constructor arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
		end
    end
    methods (Access = protected)
        
        function stepImpl(obj, y)
            x = obj.rfferx.step(y);
			y = obj.phyrx.step(x);

			% Get information for PDSCH, DM-RS and PT-RS allocations
			[pdschIndices,dmrsIndices,dmrsSymbols,...
				ptrsIndices, ptrsSym, pdschIndicesInfo] = ...
			mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
			
			for ncc=1:obj.componentCarrier

				% Perform RX beamforming by multiplying y with the RX BF vector.
				z = y(:,:,ncc)*obj.rxBF;                       

				% Demodulate the RX signal
				obj.ofdmGrid = mmwsim.nr.hOFDMDemodulate(obj.carrierConfig, z);
				
				% Estimate the Phase Noise
				obj.ofdmGrid = mmwsim.nr.PTRSEstimate(ptrsIndices, ptrsSym, obj.ofdmGrid);
						
				% Get channel estimate.
				% This is a poor channel estimate since we have not done
				% carrier and timing estimation.  But, this is OK for now.
				[chanEstGrid, obj.noiseEst, estInfo] = nrChannelEstimate(...
				obj.ofdmGrid,dmrsIndices,dmrsSymbols,...
				'CyclicPrefix',obj.carrierConfig.CyclicPrefix,...
				'CDMLengths',pdschIndicesInfo.CDMLengths);

				% Extract raw symbols and channel estimate on PDSCH
				obj.pdschSymRaw(:,ncc) = obj.ofdmGrid(pdschIndices);
				obj.pdschChanEst(:,ncc) = chanEstGrid(pdschIndices);
				obj.pdschSymEq(:,ncc) = obj.pdschSymRaw(:,ncc) ./ obj.pdschChanEst(:,ncc);
				
				
			end
        end
        
    end
end

