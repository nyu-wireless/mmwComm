classdef RFFETx < matlab.System
    % RFFETx.  Class containing the Transmit RF front-end 
    
    properties
		nonLin = false;
		phaseNoise = false;
		% Power amplifier parameters
        pa;
    end
    
    methods
        function obj = RFFETx(varargin)
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
            % power amplifier
		end
		
        function y = stepImpl(obj, x)
            % For the moment we assume ideal components at the Transmitter.
            y = x; 
		end
    end
end

