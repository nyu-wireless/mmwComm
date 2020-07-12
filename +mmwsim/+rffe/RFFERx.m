classdef RFFERx < matlab.System
    % RFFERx.  Class containing the Receive RF front-end.
	properties
		nonLin = false;
		phaseNoise = false;
		
		% LNA parameters
        lna;
        lnaMethod = 'Rapp model';     
        lnaSmooth = 1.55;	% p parameter
        lnaLinGain = 1;	% Linear gain of the LNA in dB
        lnaSat = 10;     % Limit the output signal level of the LNA
        
        % Mixer parameters
        pnoise;
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
			% Apply phase noise to waveform
			fc = 28e9;
			sr = 491.52e6;
			foffsetLog = (4.5:0.1:log10(sr/2));
			foffset = 10.^foffsetLog;
			PN_dBc_Hz = obj.PNmodelPoleZero(foffset,fc);
			obj.pnoise = comm.PhaseNoise('FrequencyOffset',foffset,'Level', PN_dBc_Hz,'SampleRate',sr);			
			
            % Create and configure a memoryless nonlinearity to model the 
            % amplifier
            obj.lna = comm.MemorylessNonlinearity;
            obj.lna.Method = obj.lnaMethod;
            obj.lna.Smoothness = obj.lnaSmooth;    
            obj.lna.LinearGain = obj.lnaLinGain;
            obj.lna.OutputSaturationLevel = obj.lnaSat;
			
			% plotNonLinearCharacteristic(obj.lna);
		end
		
        function y = stepImpl(obj, x)
			
			if obj.phaseNoise
				xpn = zeros(size(x));
				for i = 1:size(xpn, 2)
					xpn(:, i) = obj.pnoise(x(:, i));
				end
			else
				xpn = x;
			end
			
			% Apply memory-less non linearity
			if obj.nonLin
				y = zeros(size(xpn));
				for i = 1:size(xpn, 2)
					y(:, i) = obj.lna(xpn(:, i));
				end
			else
				y = xpn;
			end
		end
		
		function PN_dBC_Hz = PNmodelPoleZero(obj, f,fc)
			% Generate the phase noise characteristic in dBc/Hz for the frequency
			% offset values specified by vector f for the carrier frequency fc.
			%
			% The model used here is the multi-pole/zero model as proposed in:
			% Yinan Qi et al. "On the Phase Tracking Reference Signal (PT-RS) Design
			% for 5G New Radio (NR)". Vehicular Technology Conference (VTC-Fall), Aug
			% 2018.

			% Pole/zeros and PSD0 as specified in the paper mentioned above and in 3GPP
			% R1-163984, "Discussion on phase noise modeling". This corresponds to
			% "parameter set A" for a base frequency of 30GHz.
			fcBase = 30e9;
			fz= [1.8 2.2 40]*1e6;
			fp = [0.1 0.2 8]*1e6;
			PSD0 = -79.4;

			% Compute numerator
			num = ones(size(f));
			for ii=1:numel(fz)
				num = num.* ( 1 +(f./fz(ii)).^2);
			end

			% Compute denominator
			den = ones(size(f));
			for ii=1:numel(fz)
				den = den.* ( 1 +(f./fp(ii)).^2);
			end

			% Compute phase noise and apply a shift for carrier frequencies different
			% from the base frequency.
			PN_dBC_Hz = 10*log10(num./den) + PSD0 + 20*log10(fc/fcBase);
		end
	end
end