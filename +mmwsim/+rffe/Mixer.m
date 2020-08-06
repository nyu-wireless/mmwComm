classdef Mixer < matlab.System
	properties
		nonLinModel;	% non-linear model for the lna
		linGain = 1;	% Linear gain of the lna in dB
		satLev = 10;     % Limit the output signal level of the lna
		method = 'Rapp model';	%
		smooth =  1.55;	% p parameter		
		nonLin = true;
		phaseNoise;
		noiseFig;
		noise;
		fsamp;
	end
	
	methods
		function obj = Mixer(varargin)
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
			obj.noise = comm.ThermalNoise('NoiseMethod', 'Noise figure', ...
				'NoiseFigure', obj.noiseFig, 'SampleRate', obj.fsamp);
		end
		
		function y = stepImpl(obj, x)
			% Apply memory-less non linearity
			if obj.nonLin
				y = reshape(obj.nonLinModel(obj.noise.step(x(:))), length(x), []);
			else
				y = x;
			end
			
		end
	end			
end