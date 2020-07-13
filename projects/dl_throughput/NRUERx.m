classdef NRUERx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;  % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % ADC parameters
        nbitsADC = 0;
        adcInputVar = 1.0;
        phaseDither = true;
        
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
       
        % LNA
		satLev;
		nonLin = false;
		phaseNoise = false;
        lna;
        
		% ADC
		adc;
        
        % Carrier aggregation
        ncc = 1;
		ccFreq;
		
		% PHY
		phy;
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
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end                  
            
            % Create the LNA
            obj.lna = mmwsim.rffe.LNA('nonLin', obj.nonLin, ...
				'phaseNoise', obj.phaseNoise, 'satLev', obj.satLev);
            
            % Create the ADC
            obj.adc = mmwsim.rffe.ADC('nbits', obj.nbitsADC, 'isComplex', true, ...
			'inputVar', obj.adcInputVar, 'phaseDither', obj.phaseDither, ...
			'fsamp', obj.waveformConfig.SamplingRate);
						
            % Create the PHY
            obj.phy = mmwsim.phy.PhyRx('nbitsADC', obj.nbitsADC, 'ncc', obj.ncc, ...
				'carrierConfig', simParam.carrierConfig, ...
				'fsamp', obj.waveformConfig.SamplingRate);
        end
    end
    methods (Access = protected)
        function setupImpl(obj)
            % Set up.
            % In this case, we calibrate the ADC based on the input
            % variance
            obj.adc.optScale();
            
            % Need to add here an automatic way to calculate the center
            % frequency for each component carrier.
			
            if obj.ncc == 4
                obj.ccFreq = [-150e6, -50e6, 50e6, 150e6];
            elseif obj.ncc == 8
                obj.ccFreq = [-700e6, -500e6, -300e6, -100e6, ...
								100e6, 300e6, 500e6, 700e6];
            else
                obj.ccFreq = 0;
			end
            obj.phy.set('ccFreq', -obj.ccFreq);
        end
        
        function stepImpl(obj, y)
            
			% Pass the data from the RFFE
			y = obj.lna.step(y);

			% Pass the data from the ADC
			y = obj.adc.step(y);

			% Pass the data from the PHY
			y = obj.phy.step(y);
			
            % Loop over Component Carrier
            for cc = 1:obj.ncc
                % Get information for PDSCH, DM-RS and PT-RS allocations
                [pdschIndices,dmrsIndices,dmrsSymbols,ptrsIndices, ...
                    ptrsSym, pdschIndicesInfo] = ...
                    mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);

                % Perform RX beamforming by multiplying y with the
                % the RX BF vector.
                z = y(:,:,cc)*obj.rxBF;

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
                obj.pdschSymRaw(:, cc) = obj.ofdmGrid(pdschIndices);
                obj.pdschChanEst(:, cc) = chanEstGrid(pdschIndices);
                obj.pdschSymEq(:, cc) = obj.pdschSymRaw(:, cc) ./ obj.pdschChanEst(:, cc);
            end
        end
        
    end
end

