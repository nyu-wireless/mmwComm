classdef PDSCHSimParam < matlab.mixin.SetGet
    % PDSCHSimParam:  Parameters for the PDSCHSimParam    
    properties
                
        % Settable parameters.  These should be set in the constructor
        NLayers = 1;    % number of layers
        Modulation = 'QPSK';  
        NRB = 66;  % number of resource blocks
        SubcarrierSpacing = 120;  % SCS in kHZ
        fc = 28e9;  % carrier frequency in Hz
        
        % Carrier parameters.  Right now, we assume only one BWP
        carrierConfig; 
                 
        % PDSCH parameters
        pdschConfig;
        
        % Waveform parameters
        waveformConfig;

		% ADC resolution
		nbadc = 0;
		
		% Carrier Aggregation
		ncc;
    end
    
    methods
        function obj = PDSCHSimParam(varargin)
            % Constructor
            
            % Set parameters from constructor arguments
			if nargin >= 1
				obj.set(varargin{:});
			end

			if obj.fc == 140e9
				obj.SubcarrierSpacing = 240;
			end
			
            % Set the carrierParam settings
            obj.carrierConfig = mmwsim.nr.CarrierConfig(...
                'NRB', obj.NRB, 'SubcarrierSpacing', obj.SubcarrierSpacing, ...
                'fc', obj.fc);
            
            % Compute the waveform parameters
            % We need to convert the carrier configuration to a structure
            carrierStruct = mmwsim.nr.objtostruct( obj.carrierConfig );
            obj.waveformConfig = mmwsim.nr.hOFDMInfo(carrierStruct);
			obj.waveformConfig.SamplingRate = obj.ncc * obj.waveformConfig.SamplingRate;
			
            % PDSCH parameters.  We allocate all the RBs and all 
            % the OFDM symbols in the slot
            res.Symbols = 0;  % Reserve first symbol
            res.PRB = (0:obj.NRB-1);      % Reserve all RBs
            res.Period = 1;  % Reserve on every slot
			
            obj.pdschConfig = mmwsim.nr.PDSCHConfig(...
                'PRBSet', (0:obj.NRB-1), ...
                'SymbolSet', (0:obj.waveformConfig.SymbolsPerSlot-1), ...
                'Reserved', res, ...
				'EnablePTRS', 1,...
				'PTRSFrequencyDensity', 2, ...
				'PTRSTimeDensity', 1, ...
				'PTRSREOffset', '00'); 
		end
    end
end

