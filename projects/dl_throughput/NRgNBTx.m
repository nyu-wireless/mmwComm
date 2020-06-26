classdef NRgNBTx < matlab.System
    % 5G NR gNB transmitter class    
    properties                
        carrierConfig;  % Carrier configuration         
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grids
        ofdmGridLayer;     % Before pre-coding nsc x nsym x nlayers       
        ofdmGridAnt;       % Before pre-coding nsc x nsym x nantennas        
        
        % OFDM grid to visualize the type of symbols
        ofdmGridChan;
               
        % Transmitted data in last slots
        bits;           % TX bits
        pdschSym;       % TX symbols
        dmrsSym;        % TX data symbols
        ptrsSym;        % TX data symbols
        
        % Slot number
        Nslot = 0;
        
        % TX beamforming vector.  This is fixed.
        txBF;
        
        % Carrier aggregation
        ncc = 1;
		ccFreq;
        
        % RFFE
		nonLin = false;
		phaseNoise = false;
        rffe;
        
		% DAC
		dac;
        
        % DAC parameters
        nbitsDAC = 0;
        dacInputVar = 1.0;
        phaseDither = false;
    end
    
    properties (Constant)
        % Indices for ofdmGridChan indicating the type of symbol
        
        
    end
    
    methods
        function obj = NRgNBTx(simParam, varargin)
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
                        
            % Create the RFFE
            obj.rffe = mmwsim.rffe.RFFETx('nonLin', obj.nonLin, 'phaseNoise', obj.phaseNoise);
            
            % Create the ADC
            obj.dac = mmwsim.rffe.DAC('nbits', obj.nbitsDAC, 'isComplex', true, ...
                'inputVar', obj.dacInputVar, 'phaseDither', obj.phaseDither);
        end
    end
    methods (Access = protected)
    
        function setupImpl(obj)
            % Need to add here an automatic way to calculate the center
            % frequency for each component carrier.
            if obj.ncc == 4
                obj.ccFreq = [-150e6, -50e6, 50e6, 150e6];
            else
                obj.ccFreq = 0;
            end
            
            obj.dac.optScale();
        end
        
        function x = stepImpl(obj)
            % step implementation. Creates one slot of samples for each
            % component carrier
            
            % Set the slot number if the PDSCH config
            obj.pdschConfig.Nslot = obj.Nslot;
            
            % Loop over Component Carrier
            for cc = 1:obj.ncc

                % Create the PDSCH grid before pre-coding
                obj.ofdmGridLayer = zeros(...
                    obj.waveformConfig.NSubcarriers,...
                    obj.waveformConfig.SymbolsPerSlot, ...
                    obj.pdschConfig.NLayers);
                obj.ofdmGridChan = zeros(...
                    obj.waveformConfig.NSubcarriers,...
                    obj.waveformConfig.SymbolsPerSlot, ...
                    obj.pdschConfig.NLayers);

                % Get information for PDSCH and DM-RS allocations
                [pdschIndices, dmrsIndices, obj.dmrsSym(:, cc), ...
                    ptrsIndices, obj.ptrsSym(:, cc), pdschIndicesInfo] = ...
                    mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);

                % Generate random bits
                bitsPerSym = mmwsim.nr.NRConst.bitsPerSym(obj.pdschConfig.Modulation);
                nsym = length(pdschIndices);
                nbits = bitsPerSym * nsym;
                obj.bits = randi([0 1], nbits, 1);

                % Modulate the bits to symbols
                M = mmwsim.nr.NRConst.modOrder(obj.pdschConfig.Modulation);
                obj.pdschSym(:, cc) = qammod(obj.bits,M,'InputType','bit',...
                    'UnitAveragePower',true);

                % Map symbols to OFDM grid
                obj.ofdmGridLayer(pdschIndices) = obj.pdschSym(:, cc);
                obj.ofdmGridLayer(dmrsIndices) = obj.dmrsSym(:, cc);
                obj.ofdmGridLayer(ptrsIndices) = obj.ptrsSym(:, cc);


                % Fill the channel with labels of the channels.
                % This is just for visualization
                obj.ofdmGridChan(pdschIndices) = 1;
                obj.ofdmGridChan(dmrsIndices) = 2; 
                obj.ofdmGridChan(ptrsIndices) = 3;   

                % Perform the OFDM modulation
                xlayer = mmwsim.nr.hOFDMModulate(obj.carrierConfig, obj.ofdmGridLayer);

                % Perform the TX beamforming.  At this point, 
                % xlayer will be an nsamp x 1 vector.  Use the TX beamforming
                % vector, obj.txBF, to map this to a  nsamp x nant matrix
                % where nant is the number of TX antennas.
                xbf = xlayer*obj.txBF.';
                
                % Upsample the Tx signal to the total bandwidth
                xup = resample(xbf, obj.ncc, 1)/obj.ncc;
                
                % Use a NCO to shift the component carrier
				xnco = mmwsim.nr.hCarrierAggregationModulate(xup, obj.dac.fsamp, obj.ccFreq(cc));
                
                % Add the 
				x(:,:,cc) = xnco;
            end
            x = sum(x, 3);
            % Increment the slot number
            obj.Nslot = obj.Nslot + 1;
        end        
      
    end
end

