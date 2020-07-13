classdef PhyRx < matlab.System
	% PhyRx. The physical layer processing after the A/D
	
	properties
        carrierConfig;  % Carrier configuration
		nbitsADC = 0;	% 0: Infinite processing
						% >0: Fixed-point processing
		
		% CC Filter design specifications
		Fp;				% pass-band frequency
		Fst;			% stop-band frequency
		Ap = 1;			% pass-band ripple (dB)
		Ast = 50;		% stop-band attenuation (dB)
		nbcoeff = 6;	% number of bits for the filter coefficient
		nfilt;			% filter order
		bfilt;			% num coeffs
		afilt;			% denom coeffs
		hd;				% filter design
		
        % Carrier aggregation
        ncc = 1;
		ccFreq;
		
		% ADC
		fsamp;
	end
	
	methods
		function obj = PhyRx(varargin)
            % Constructor
            
            % Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
		end
	end
	
	methods (Access = protected)
        function setupImpl(obj)
			if obj.nbitsADC > 0
				% Find the effective signal bandwidth:
				%
				% * ResourceBlocks * 12 * SubCarrierSpacing
				fsig = obj.carrierConfig.NRB * 12 * ...
					obj.carrierConfig.SubcarrierSpacing * 1e3;
				obj.Fp = fsig/obj.fsamp;
				obj.Fst = 1/obj.ncc;

				% Design a Fixed-Point Filter
				spec = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
					obj.Fp, obj.Fst, obj.Ap, obj.Ast);
				f = design(spec, 'minphase', false, 'SystemObject', true);

				% measure(f)
				coef = coeffs(f);

				bq = fi(coef.Numerator, true, obj.nbcoeff, ...
					'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
				L = bq.FractionLength;
				bsc = coef.Numerator*2^L;

				obj.hd = dfilt.dffir(bsc);

				obj.hd.Arithmetic = 'fixed';
				obj.hd.CoeffWordLength = obj.nbcoeff;
				all(obj.hd.Numerator == round(obj.hd.Numerator));

				% Integer real input from ADC with nbitsADC resolution
				obj.hd.InputWordLength = obj.nbitsADC;
				obj.hd.InputFracLength = 0;

				obj.bfilt = obj.hd.Numerator;
				obj.afilt = 1;
			end
		end
		
        function y = stepImpl(obj, x)
			if obj.ncc == 1
				% We have a single component carrier.
				y = x;
			else
				y = zeros(size(x,1)/obj.ncc, size(x,2), obj.ncc);
				for cc = 1:obj.ncc
					% Use a numerically controlled oscillator to bring the 
					% component carrier of interest to the middle of the
					% spectrum
					xnco = mmwsim.nr.hCarrierAggregationModulate(x, obj.fsamp, obj.ccFreq(cc));

					% filter the input signal
					if obj.nbitsADC > 0
						xfilt = filter(obj.bfilt, obj.afilt, xnco);
					else
						xfilt = lowpass(xnco, 0.1934,'ImpulseResponse','fir','Steepness',0.95);
					end
					
					% Downsample the filted signal
					y(:,:,cc) = resample(xfilt, 1, obj.ncc);
				end
			end
		end
	end
end