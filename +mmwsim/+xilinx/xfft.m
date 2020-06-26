classdef xfft < matlab.System
	% FFT. Bit accurate FFT model from Xilinx
	
	properties
		generics
	end
	
	methods
		function obj = xfft(varargin)
			% Constructor

			% Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end

			% Compile the fft model
		end
	end
	
	methods (Access = protected)
		function y = stepImpl(obj, x)
%			y = xfft_v9_0_bitacc_mex(obj.generics, 3, quantize(rx_fft.q,rx_fir_op_core(:,n)), 0, 0);
		end
	end
end

