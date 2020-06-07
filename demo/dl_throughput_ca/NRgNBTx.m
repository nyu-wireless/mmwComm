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
        ptrsSym;		% TX phase tracking symbols
				
        % Slot number
        Nslot = 0;
        
        % TX beamforming vector.  This is fixed.
        txBF;
		
		% DAC calibration parameters
		nscal = 1e4;		% number of samples used for calibration
		nbdac = 4;			% number of bits in the DAC
		mseQ;				% quantizer MSE
		aq;					% optimizer quantizer scale value
		autoscale = true;	% scale input (as if there is an oracle AGC)
		fsamp = 491.52e6	% sample frequency
		
		% Component Carrier Aggregation
		componentCarrier = 4;
		fc = [-150e6, -50e6, 50e6, 150e6];
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
        end
    end
    methods (Access = protected)
    
        function x = stepImpl(obj)
            % step implementation.  Creates one slot of samples
            x = zeros(61632, 64);
			
			for ncc=1:4
				% Set the slot number if the PDSCH config
				obj.pdschConfig.Nslot = obj.Nslot;

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
				[pdschIndices,dmrsIndices,obj.dmrsSym(:,ncc),...
					ptrsIndices,obj.ptrsSym(:,ncc), pdschIndicesInfo] = ...
					mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);

				% Generate random bits
				bitsPerSym = mmwsim.nr.NRConst.bitsPerSym(obj.pdschConfig.Modulation);
				nsym = length(pdschIndices);
				nbits = bitsPerSym * nsym;
				obj.bits = randi([0 1], nbits, 1);

				% Modulate the bits to symbols
				M = mmwsim.nr.NRConst.modOrder(obj.pdschConfig.Modulation);
				obj.pdschSym(:, ncc) = qammod(obj.bits,M,'InputType','bit',...
					'UnitAveragePower',true);

				% Map symbols to OFDM grid
				obj.ofdmGridLayer(pdschIndices) = obj.pdschSym(:, ncc);
				obj.ofdmGridLayer(dmrsIndices) = obj.dmrsSym(:, ncc);
				obj.ofdmGridLayer(ptrsIndices) = obj.ptrsSym(:, ncc);

				% Fill the channel with labels of the channels.
				% This is just for visualization
				obj.ofdmGridChan(pdschIndices) = 1;
				obj.ofdmGridChan(dmrsIndices) = 2;   
				obj.ofdmGridChan(ptrsIndices) = 3;    

				% Perform the OFDM modulation
				xlayer = mmwsim.nr.hOFDMModulate(obj.carrierConfig, obj.ofdmGridLayer);

				% Perform the TX beamforming.  At this point, xlayer will be 
				% an nsamp x 1 vector.  Use the TX beamforming vector, 
				% obj.txBF, to map this to a  nsamp x nant matrix where nant is
				% the number of TX antennas.
				xbf = xlayer*obj.txBF.';
				
				xup = resample(xbf, obj.componentCarrier, 1)/obj.componentCarrier;
				xnco = mmwsim.nr.hCarrierAggregationModulate(xup, obj.fsamp, obj.fc(ncc));
				x = x + xnco;
			end
			
            % Increment the slot number
            obj.Nslot = obj.Nslot + 1;
		end
      
    end
end

