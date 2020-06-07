classdef RFFERx < matlab.System
    % RFFERx.  Class containing the RF front-end of the UE.
	properties		        
		% ADC calibration parameters
		nscal = 1e4;		% number of samples used for calibration
		nbadc = 6;			% number of bits in the ADC
		mseQ;				% quantizer MSE
		aq;					% optimizer quantizer scale value
		autoscale = true;	% scale input (as if there is an oracle AGC)
		fsamp = 491.52e6;	% sample frequency in Hz
	end
	
	methods
		function obj = RFFERx(varargin)
            % Constructor
			
            % Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
		end
	end
		
    methods (Access = protected)
        function setupImpl(obj)
			% Optimize quantizer scale value
			% Find the optimal scale level for the quantizer assuming
			% the input has variance 1
			if obj.nbadc > 0
				[obj.aq,obj.mseQ] = mmwsim.rffe.Quant.optScale(obj.nbadc, obj.nscal);
				obj.mseQ = 10*log10(obj.mseQ);                          
			end
		end
		
        function y = stepImpl(obj, x)
			% Autoscale.  This would be done by the AGC.
			% Note that in the complex case, the input is scaled
			% unit variance per I/Q.
			if obj.autoscale
				xvar = mean(abs(x).^2);
				if isreal(x)
					x = x ./ sqrt(xvar);
				else
					x = x ./ sqrt(xvar/2);
				end
			end
			
			% Quantize the input.
			xq = mmwsim.rffe.Quant.qsat(x, obj.nbadc, obj.aq);
			if obj.nbadc > 0
				bq = fi(xq, true, obj.nbadc, 'FractionLength', 0, ...
					'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
				y = double(bq);
			else
				y = xq;
			end
		end
	end
end