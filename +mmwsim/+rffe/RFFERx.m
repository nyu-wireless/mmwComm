classdef RFFERx < matlab.System
    % RFFERx.  Class containing the Receive RF front-end.
	properties
		isLinear = false;	% include the rffe non-linearities
		
		% Simulation parameters
		nrx = 16;			% number of RX antennas
		noiseTemp = 290;	% noise temperature in K
		ktVar = 1.0;		% thermal noise per sample
		EkT;				% thermal noise energy
		xvar = 1;			% variance of the tx symbols
			
		lnaGain;			% LNA gain in dB
		lnaNF;				% LNA noise figure in dB
		lnaAmpLut;			% LNA fund tone power curve
		lnaPower;			% LNA power in mW
		
		mixGain;			% mixer gain in dB
		mixNF;				% mixer noise figure in dB
		mixAmpLut;			% mixer fund tone power curve
		mixPLO;				% local oscillator power
		mixPower;			% mixer power in mW
		
		nbits = 4;			% ADC bits per dimension
		adcFOM = 65e-12;	% FOM of the ADC 65 fJ/conv
		fsamp = 1e9;		% sample rate
		
		elem;				% rffe elements
	end
	
	methods
		function obj = RFFERx(varargin)
            % Constructor
			
            % Set key-value pair arguments
			if nargin >= 1
				obj.set(varargin{:});
			end
			
			obj.EkT = physconst('Boltzman')*obj.noiseTemp;
		end
		
		function NF = nf(obj)
			% Calculate the effective noise figure
			NF = 10*log10(10^(0.1*obj.lnaNF) + 10^(-0.1*obj.lnaGain)*(10^(0.1*obj.mixNF)-1));
		end
		
		function P = power(obj)
			% Calculate the power consumption in mW
			if obj.nbits == 0
				% For inf-bit ADCs use 12-bit for power calculation.
				nb = 12;
			else
				nb = obj.nbits;
			end
			
			ndriverTest = 2.^(0:log2(obj.nrx));	% number of LO drivers
			ndriver = length(ndriverTest);
			
			% Initialize the power consumption vector.
			P = zeros(ndriver, 1);
			
			% In a fully-digital receiver we have an LNA and a Mixer per
			% antenna.
			P = P + obj.nrx*obj.lnaPower;
			P = P + obj.nrx*obj.mixPower;
			
			% Find the LO power for each configuration
			L3dB = 0.5;					% loss in dB
			etaOpt = 0.25;				% maximum power-added efficiency (PAE)
			Popt = 10^(0.1*10);			% maximum output power in mW
			Pin = 10^(0.1*(-5));		% input power of the LO signal to the power divider in mW
			Pmul = 43;					% power consumption of the LO generation network in mW
			PLO = 10^(0.1*obj.mixPLO);	% power required at the input of each mixer in mW
			
			loPower = zeros(ndriver,1);
			
			for idriver = 1:ndriver
				Nd = ndriverTest(idriver);
				Pout = PLO*(obj.nrx/Nd)*10^(0.1*L3dB*log2(obj.nrx/Nd));
				Pamp = max(0, Pout - Pin);
				if Pamp > 0
					eta = etaOpt*2/((Pamp/Popt)+(Popt/Pamp));
					Pdriver = Pamp/eta;
				else
					Pdriver = 0;
				end
				loPower(idriver) = Nd*(Pdriver+Pmul);
			end
			P = P + loPower;
			
			% ADC Power
			P = P + obj.nrx*2*obj.fsamp*(2^nb)*obj.adcFOM;
		end
	end
		
    methods (Access = protected)
        function setupImpl(obj)
			
			% Thermal noise at the antenna.
			% Note:  Due to the convention in the ThermalNoise object in 
			% MATLAB, you must set SampleRate=1 to get the correct noise 
			% energy per sample of kT.
			tn0 = comm.ThermalNoise('SampleRate', 1, 'NoiseTemperature', ...
				obj.noiseTemp);
			tn = MultiInput(tn0, obj.nrx);
				
			% LNA thermal noise
			lnaNoise0 = comm.ThermalNoise('NoiseMethod', 'Noise figure', ...
				'NoiseFigure', obj.lnaNF, 'SampleRate', 1);
			lnaNoise = MultiInput(lnaNoise0, obj.nrx);
			
			% LNA gain and nonlinearity
			lnaAmp0 = comm.MemorylessNonlinearity('Method', 'Lookup table', ...
				'Table', obj.lnaAmpLut);
			lnaAmp = MultiInput(lnaAmp0, obj.nrx);

			% Mixer thermal noise
			mixNoise0 = comm.ThermalNoise('NoiseMethod', 'Noise figure', ...
				'NoiseFigure', obj.mixNF, 'SampleRate', 1);
			mixNoise = MultiInput(mixNoise0, obj.nrx);
			
			% Mixer gain and nonlinearity
			mixAmp0 = comm.MemorylessNonlinearity('Method', 'Lookup table', ...
				'Table', obj.mixAmpLut);
			mixAmp = MultiInput(mixAmp0, obj.nrx);

			% Baseband AGC used to adjust the input level to the ADC-
			% This would be performed via a controllable baseband amplifier
			bbAGC = mmwsim.rffe.AutoScale('meth', 'MatchTgt');
						
			% ADC
			adc = mmwsim.rffe.ADC('nbits', obj.nbits);
			
			% Find optimal input target for the ADC
			% Full scale value
			adcFS = max(adc.stepSize*(2^(obj.nbits-1)-1), 1);
			EsFS = 2*adcFS^2;

			% Test the values at some backoff from full scale
			bkfTest = linspace(-30,0,100)';
			EsTest = EsFS*10.^(0.1*bkfTest);

			% Compute the SNR 
			snr = bbAGC.compSnrEs(EsTest, adc);

			% Select the input level with the maximum SNR
			[~, im] = max(snr);
			bbAGC.set('EsTgt', EsTest(im));

			obj.elem = {tn, lnaNoise, lnaAmp, mixNoise, mixAmp, bbAGC, adc};
		end
		
        function y = stepImpl(obj, x)	
			% Find the number of RFFE elements
			nstage = length(obj.elem);
			
            % Add thermal noise
            x = obj.elem{1}.step(x);
            
			if ~obj.isLinear
				% Apply memory-less non linearity 
				for i = 2:nstage-2
					x = obj.elem{i}.step(x);
				end
			end
			
			% The last stage of the RFFE is the ADC
			x = obj.elem{nstage-1}.step(x);
            
			% The last stage of the RFFE is the ADC
			y = obj.elem{nstage}.step(x);
		end
	end
end