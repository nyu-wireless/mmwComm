classdef xfir < matlab.System
	% FIR. Bit accurate FIR model from Xilinx
	
	properties
		generics
	end
	
	methods
		function obj = xfir(varargin)
			% Constructor

			% Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end

			% Compile the fir model
		end
	end
	
	methods (Access = protected)
		function y = stepImpl(obj, x)
%			y = xfir_v7_2_bitacc_mex(...);
		end
	end
end

