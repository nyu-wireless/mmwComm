classdef NRUERx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;  % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grid b
        ofdmGrid;     % Before pre-coding nsc x nsym x nlayers
        
        % Channel and noise estimate
        noiseEst;
        
        % Recived symbols
        pdschChanEst;   % Channel estimate on the PDSCH
        pdschSymRaw;    % Raw symbols before equalization
        pdschSymEq;     % Equalized symbols
        dmrsSym;        % DM-RS Reference symbols
        
        % Slot number
        Nslot = 0;
        
        % RX beamforming vector.  
        rxBF;
        
        % Timing offset
        offset;
        
		% ADC calibration parameters
		nscal = 1e4;		% number of samples used for calibration
		nbadc = 4;			% number of bits in the ADC
		mseQ;				% quantizer MSE
		aq;					% optimizer quantizer scale value
		autoscale = true;	% scale input (as if there is an oracle AGC)
		fsamp = 491.52e6	% sample frequency
		
		% Filter design specifications
		Fp;				% pass-band frequency
		Fst;			% stop-band frequency
		Ap = 0.1;			% pass-band ripple (dB)
		Ast = 50;		% stop-band attenuation (dB)
		nbcoeff = 6;	% number of bits for the filter coefficient
		nfilt;			% filter order
		bfilt;			% num coeffs
		afilt;			% denom coeffs
		hd;				% filter design
		
		% Component Carrier Aggregation
		componentCarrier = 4;
		fc = [150e6, 50e6, -50e6, -150e6];
		firFilter;
    end
    
    methods
        function obj = NRUERx(simParam, varargin)
            % Constructor
           
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the 
            % objects need to be converted to older structures.
            obj.carrierConfig = mmwsim.nr.objtostruct( simParam.carrierConfig );
            obj.pdschConfig = mmwsim.nr.objtostruct( simParam.pdschConfig );
            obj.waveformConfig = mmwsim.nr.objtostruct( simParam.waveformConfig );
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
			end                        
            
			% Design the downsampling component-carrier filter
			obj.designFilt()
		end
		
		function designFilt(obj)
			fsig = obj.carrierConfig.NRB*12*obj.carrierConfig.SubcarrierSpacing*1e3;
			obj.Fp = fsig/obj.fsamp;
			obj.Fst = 1/obj.componentCarrier;
			
% 			spec = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
% 				obj.Fp, obj.Fst, obj.Ap, obj.Ast);
% 			f = design(spec, 'minphase', false, 'SystemObject', true);
% 			coef = coeffs(f);
			Fpass = 0.1934;           % Passband Frequency
			Fstop = 0.25;             % Stopband Frequency
			Dpass = 0.057501127785;   % Passband Ripple
			Dstop = 0.0031622776602;  % Stopband Attenuation
			flag  = 'scale';          % Sampling Flag

			% Calculate the order from the parameters using KAISERORD.
			[N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop], [1 0], [Dstop Dpass]);

			% Calculate the coefficients using the FIR1 function.
			b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
			obj.bfilt = b;
			obj.afilt = 1;

			% The response of the designed filter is displayed.
			% fvtool(obj.bfilt, obj.afilt,'Analysis','freq');
		end
    end
    methods (Access = protected)
        
        function stepImpl(obj, y)
			% Get information for PDSCH, DM-RS and PT-RS allocations
			[pdschIndices,dmrsIndices,dmrsSymbols,...
				ptrsIndices, ptrsSym, pdschIndicesInfo] = ...
			mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
			
			for ncc=1:obj.componentCarrier
				% Use a numerically controlled oscillator to bring the 
				% component carrier of interest to the middle of the
				% spectrum
				ynco = mmwsim.nr.hCarrierAggregationModulate(y, obj.fsamp, obj.fc(ncc));
				
				% filter the input signal
				yfilt = lowpass(ynco, 0.1934,'ImpulseResponse','fir','Steepness',0.95);
				
				% Downsample the filted signal
				ydown = resample(yfilt, 1, obj.componentCarrier);

				% Perform RX beamforming by multiplying y with the RX BF vector.
				z = ydown*obj.rxBF;                       

				% Demodulate the RX signal
				obj.ofdmGrid = mmwsim.nr.hOFDMDemodulate(obj.carrierConfig, z);
				
				% Estimate the Phase Noise
				obj.ofdmGrid = mmwsim.nr.PTRSEstimate(ptrsIndices, ptrsSym, obj.ofdmGrid);
						
				% Get channel estimate.
				% This is a poor channel estimate since we have not done
				% carrier and timing estimation.  But, this is OK for now.
				[chanEstGrid, obj.noiseEst, estInfo] = nrChannelEstimate(...
				obj.ofdmGrid,dmrsIndices,dmrsSymbols,...
				'CyclicPrefix',obj.carrierConfig.CyclicPrefix,...
				'CDMLengths',pdschIndicesInfo.CDMLengths);

				% Extract raw symbols and channel estimate on PDSCH
				obj.pdschSymRaw(:,ncc) = obj.ofdmGrid(pdschIndices);
				obj.pdschChanEst(:,ncc) = chanEstGrid(pdschIndices);
				obj.pdschSymEq(:,ncc) = obj.pdschSymRaw(:,ncc) ./ obj.pdschChanEst(:,ncc);
				
				
			end
        end
        
    end
end

