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
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end                        
            
        end
    end
    methods (Access = protected)
        
        function stepImpl(obj, y)
            y = obj.rfferx.step(y);
			
            % Get information for PDSCH, DM-RS and PT-RS allocations
            [pdschIndices,dmrsIndices,dmrsSymbols,ptrsIndices, ...
				ptrsSym, pdschIndicesInfo] = ...
                mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
			
            % Perform RX beamforming by multiplying y with the
            % the RX BF vector.
            z = y*obj.rxBF;
            
            % Demodulate the RX signal
            obj.ofdmGrid = mmwsim.nr.hOFDMDemodulate(obj.carrierConfig, z);
			
			% Estimate the Phase Noise
			obj.ofdmGrid = mmwsim.nr.PTRSEstimate(ptrsIndices, ptrsSym, obj.ofdmGrid);
			
            % Get channel estimate.
            % This is a poor channel estimate since we have not done
            % carrier and timing estimation.  But, this is OK for now.
            [chanEstGrid, obj.noiseEst] = nrChannelEstimate(...
                obj.ofdmGrid,dmrsIndices,dmrsSymbols,...
                'CyclicPrefix',obj.carrierConfig.CyclicPrefix,...
                'CDMLengths',pdschIndicesInfo.CDMLengths);
			
            % Extract raw symbols and channel estimate on PDSCH
            obj.pdschSymRaw = obj.ofdmGrid(pdschIndices);
            obj.pdschChanEst = chanEstGrid(pdschIndices);
            obj.pdschSymEq = obj.pdschSymRaw ./ obj.pdschChanEst;
            
        end
        
    end
end

