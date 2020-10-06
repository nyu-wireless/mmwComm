classdef AutoScale < matlab.System
	% Simple scaling to model an AGC

	properties
		% Energy per sample target
		EsTgt = 1;

		% Number of inputs for calibration
		ncal = 1e4;

		% Scale method
		% 'None':  No scaling
		% 'MatchTgt':  Scaling to match target
		% 'AttnOnly':  Attenuate only
		methVals = {'None', 'MatchTgt', 'AttnOnly'};
		meth = 'MatchTgt';

		% Last scale level in dB
		scaledB;
	end
    
	% Quantizer methods.  All methods are static
	methods
		
		function obj = AutoScale(varargin)
			% Constructor
			% Set parameters from constructor arguments
			if nargin >= 1
				obj.set(varargin{:});
			end

			% Check arguments
			nmeth = length(obj.methVals);
			found = false;
			
			for i = 1:nmeth
				if strcmp(obj.meth, obj.methVals{i})
					found = true;
				end
			end
			
			if ~found
				error('Method %s not valid', obj.meth);
			end          
		end
		
		function snr = compSnrEs(obj, EsTest, sys)
			% Finds the distortion of a system, sys, for different
			% input levels, EsTest.  For each input level, it outputs
			% the MSE in dB.  This can be used to find the optimal
			% scaling input.  The object sys must have a method
			% y = sys.step(x) that can operate on a vector x.

			% Generate unit input test
			x = (randn(obj.ncal,1) + 1i*randn(obj.ncal,1))*sqrt(1/2);
			ntest = length(EsTest);
			snr = zeros(ntest,1);
			
			for i = 1:ntest
				% Scale input
				xs = x*sqrt(EsTest(i));
				ys = sys.step(xs);

				% Find linear model ys = a*xs + d,  d ~ CN(0,dvar)
				a = mean(conj(ys).*xs)/mean(abs(xs).^2);
				dvar = mean(abs(ys - a*xs).^2);

				% Measure the SNR
				snr(i) = 10*log10(abs(a).^2*mean(abs(xs).^2) / dvar);
			end
		end
	end

	methods (Access = protected)
		function xs = stepImpl(obj, x)
			% Step function:  Performs the scaling
			Ex = mean(abs(x).^2,'all');
			xs = x;
			if Ex == 0                
				return
			end

			% Compute the scale
			scale = sqrt(obj.EsTgt/Ex);

			% Perform the scale            
			if strcmp(obj.meth, 'MatchTgt')
				xs = scale*x;
			elseif strcmp(obj.meth, 'AttnOnly') && (scale < 1)
				xs = scale*x;          
			end

			% Record the scale in dB
			obj.scaledB = 10*log10(scale);
		end
	end
end