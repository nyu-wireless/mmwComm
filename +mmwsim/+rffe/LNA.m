classdef LNA < matlab.System
	properties
		nonLinModel;	% non-linear model for the lna
		linGain = 1;	% Linear gain of the lna in dB
		satLev = 10;     % Limit the output signal level of the lna
		method = 'Rapp model';	%
		smooth =  1.55;	% p parameter		
		nonLin = true;
		phaseNoise;
	end
	
	methods
		function obj = LNA(varargin)
            % Constructor
			
            % Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
		end
	end	
	
    methods (Access = protected)
        function setupImpl(obj)
            % Create and configure a memoryless nonlinearity to model the 
            % amplifier
            obj.nonLinModel = comm.MemorylessNonlinearity;
            obj.nonLinModel.Method = obj.method;
            obj.nonLinModel.Smoothness = obj.smooth;    
            obj.nonLinModel.LinearGain = obj.linGain;
            obj.nonLinModel.OutputSaturationLevel = obj.satLev;
		end
		
		function y = stepImpl(obj, x)
			% Apply memory-less non linearity
			if obj.nonLin
				y = reshape(obj.nonLinModel(x(:)), length(x), []);
			else
				y = x;
			end
		end
	end			
end

