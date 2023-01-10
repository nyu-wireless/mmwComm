classdef RFFERx < matlab.System
    % RFFERx.  Class containing the Receive RF front-end.
    properties
        isLinear = false;	% include the rffe non-linearities
        
        % Simulation parameters
        nrx = 16;         % number of RX antennas
        noiseTemp = 290;  % noise temperature in K
        ktVar = 1.0;      % thermal noise per sample
        EkT;              % thermal noise energy
        xvar = 1;         % variance of the tx symbols
        
        lna;
        mixer;
        ps;
        nstages = 1;
        
        nbits = 4;			% ADC bits per dimension
        adcFOM = 65e-12;	% FOM of the ADC 65 fJ/conv
        SampleRate;		% sample rate
        
        elem;				% rffe elements
        isFD = true;
        psLoss = 16;        % phase shifter loss in dB
        iir;
        isFilter = false;
    end
    
    methods
        function obj = RFFERx(varargin)
            % Constructor
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            obj.EkT = physconst('Boltzman')*obj.noiseTemp;
            
            obj.iir = mmwsim.rffe.IIRFilter();
        end
        
        function NF = nf(obj)
            % Calculate the effective noise figure
            if obj.isFD
                % LNA
                F0 = 10^(0.1*obj.lna.NF);
                G0 = 10^(0.1*obj.lna.Gain);
                
                % Mixer
                F1 = 10^(0.1*obj.mixer.NF);
                G1 = 10^(0.1*obj.mixer.Gain);
                
                F = F0 + (F1-1)/G0;
                
                NF = 10*log10(F);
            elseif obj.nstages==1
                % LNA 0
                F0 = 10^(0.1*obj.lna.NF);
                G0 = 10^(0.1*obj.lna.Gain);
                
                % LNA 1
                F1 = 1;
                G1 = 1;
                
                % Phase shihfter
                F2 = 10^(0.1*obj.ps.NF);
                G2 = 10^(0.1*obj.ps.Gain);
                
                % Combiner
                F3 = 10.^(0.1*log2(obj.nrx));
                G3 = 1/F3;
                
                % Mixer
                F4 = 10^(0.1*obj.mixer.NF);
                G4 = 10^(0.1*obj.mixer.Gain);
                
                F = F0 + (F1 - 1)/G0 + (F2 - 1)/(G0 * G1) + ...
                    (F3 - 1)/(G0 * G1 * G2) + (F4 - 1)/(G0 * G1 * G2 * G3);
                NF = 10*log10(F);
            elseif obj.nstages==2
                % LNA 0
                F0 = 10^(0.1*obj.lna.NF);
                G0 = 10^(0.1*obj.lna.Gain);
                
                % LNA 1
                F1 = 10^(0.1*obj.lna.NF);
                G1 = 10^(0.1*obj.lna.Gain);
                
                % Phase shihfter
                F2 = 10^(0.1*obj.ps.NF);
                G2 = 10^(0.1*obj.ps.Gain);
                
                % Combiner
                F3 = 10.^(0.1*log2(obj.nrx));
                G3 = 1/F3;
                
                % Mixer
                F4 = 10^(0.1*obj.mixer.NF);
                G4 = 10^(0.1*obj.mixer.Gain);
                
                F = F0 + (F1 - 1)/G0 + (F2 - 1)/(G0 * G1) + ...
                    (F3 - 1)/(G0 * G1 * G2) + (F4 - 1)/(G0 * G1 * G2 * G3);
                NF = 10*log10(F);
            elseif obj.nstages==3
                % LNA 0
                F0 = 10^(0.1*obj.lna.NF);
                G0 = 10^(0.1*obj.lna.Gain);
                
                % LNA 1
                F1 = 10^(0.1*obj.lna.NF);
                G1 = 10^(0.1*obj.lna.Gain);
                % LNA 2
                F2 = 10^(0.1*obj.lna.NF);
                G2 = 10^(0.1*obj.lna.Gain);
                
                % Phase shihfter
                F3 = 10^(0.1*obj.ps.NF);
                G3 = 10^(0.1*obj.ps.Gain);
                
                % Combiner
                F4 = 10.^(0.1*log2(obj.nrx));
                G4 = 1/F4;
                
                % Mixer
                F5 = 10^(0.1*obj.mixer.NF);
                G5 = 10^(0.1*obj.mixer.Gain);
                
                F = F0 + (F1 - 1)/G0 + (F2 - 1)/(G0 * G1) + ...
                    (F3 - 1)/(G0 * G1 * G2) + ...
                    (F4 - 1)/(G0 * G1 * G2 * G3) + ...
                    (F5 - 1)/(G0 * G1 * G2 * G3 * G4);
                NF = 10*log10(F);
            end
        end
        
        function P = power(obj)
            % Calculate the power consumption in mW
            if obj.nbits == 0
                % For inf-bit ADCs use 12-bit for power calculation.
                nb = 12;
            else
                nb = obj.nbits;
            end
            
            Pmul = 43;  % power consumption of the LO generation network in mW
            
            if obj.isFD
                ndriverTest = 2.^(0:log2(obj.nrx));	% number of LO drivers
                ndriver = length(ndriverTest);

                % Initialize the power consumption vector.
                P = zeros(ndriver, 1);

                % In a fully-digital receiver we have an LNA and a Mixer per
                % antenna.
                P = P + obj.nstages*obj.nrx*obj.lna.Power;
                P = P + obj.nrx*obj.mixer.Power;

                % Find the LO power for each configuration
                L3dB = 0.5;					% loss in dB
                etaOpt = 0.25;				% maximum power-added efficiency (PAE)
                Popt = 10^(0.1*10);			% maximum output power in mW
                Pin = 10^(0.1*(-5));		% input power of the LO signal to the power divider in mW
                
                PLO = 10^(0.1*obj.mixer.PLO);	% power required at the input of each mixer in mW

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
                P = P + obj.nrx*2*obj.SampleRate*(2^nb)*obj.adcFOM;
            else
                P = obj.nstages * obj.nrx*obj.lna.Power;          % LNA
                P = P + obj.nrx*obj.ps.Power;           % Phase shifter
                P = P + obj.mixer.Power;                % Mixer
                P = P + 10^(0.1*obj.mixer.PLO) + Pmul;  % LO
                P = P + obj.nrx*2*obj.SampleRate*(2^nb)*obj.adcFOM; % ADC
            end
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            % Thermal noise at the antenna
            tn0 = comm.ThermalNoise(...
                'SampleRate', obj.SampleRate, ...
                'NoiseTemperature', obj.noiseTemp);
            tn = MultiInput(tn0, obj.nrx);
            
            % LNA thermal noise
            lnaNoise0 = comm.ThermalNoise(...
                'SampleRate', obj.SampleRate, ...
                'NoiseMethod', 'Noise figure', ...
                'NoiseFigure', obj.lna.NF);
            lnaNoise = MultiInput(lnaNoise0, obj.nrx);
            
            % LNA gain and nonlinearity
            lnaAmp0 = comm.MemorylessNonlinearity(...
                'Method', 'Lookup table', ...
                'Table', obj.lna.AmpLut);
            lnaAmp = MultiInput(lnaAmp0, obj.nrx);
            
            if obj.isFD
                psNoise = lnaNoise;
                psAmp = lnaAmp;
            else
                % Phase shifter thermal noise
                psNoise0 = comm.ThermalNoise(...
                    'SampleRate', obj.SampleRate, ...
                    'NoiseMethod', 'Noise figure', ...
                    'NoiseFigure', obj.ps.NF);
                psNoise = MultiInput(psNoise0, obj.nrx);

                % Phase shifter gain and nonlinearity
                psAmp0 = comm.MemorylessNonlinearity(...
                    'Method', 'Lookup table', ...
                    'Table', obj.ps.AmpLut);
                psAmp = MultiInput(psAmp0, obj.nrx);
            end
            
            combNoise0 = comm.ThermalNoise(...
                'SampleRate', obj.SampleRate, ...
                'NoiseMethod', 'Noise figure', ...
                'NoiseFigure', log2(obj.nrx));
            combNoise = MultiInput(combNoise0, obj.nrx);
            
            % Mixer thermal noise
            mixNoise0 = comm.ThermalNoise(......
                'SampleRate', obj.SampleRate, ...
                'NoiseMethod', 'Noise figure', ...
                'NoiseFigure', obj.mixer.NF);
            
            if obj.isFD
                mixNoise = MultiInput(mixNoise0, obj.nrx);
            else
                mixNoise = MultiInput(mixNoise0, 1);
            end
            
            % Mixer gain and nonlinearity
            mixAmp0 = comm.MemorylessNonlinearity(...
                'Method', 'Lookup table', ...
                'Table', obj.mixer.AmpLut);
            if obj.isFD
                mixAmp = MultiInput(mixAmp0, obj.nrx);
            else
                mixAmp = MultiInput(mixAmp0, 1);
            end
            
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
            
            obj.elem = {tn, lnaNoise, lnaAmp, psNoise, psAmp, ...
                combNoise, mixNoise, mixAmp, bbAGC, adc};
        end
        
        function y = stepImpl(obj, x, wrx)
            % Find the number of RFFE elements
            nstage = length(obj.elem);
            
            % Add thermal noise, 
            x = obj.elem{1}.step(x);
            
            if ~obj.isLinear
                for i=1:obj.nstages
                    x = obj.elem{2}.step(x);  % LNA thermal noise
                    x = obj.elem{3}.step(x);  % LNA AM/AM and AM/PM
                end
            end
            
            % if the receiver uses analog or hybrid beamforming apply the
            % phase shifter loss.
            if ~obj.isFD
                if ~obj.isLinear
                    x = obj.elem{4}.step(x);  % Phase shifter thermal noise
                    x = obj.elem{5}.step(x);  % Phase shifter AMAM and AMPM
                end
                % Combiner
                x = obj.elem{6}.step(x); % Combiner noise
                x = sqrt(1/db2pow(log2(obj.nrx))) * x; % Combiner loss
                
                x = x * wrx;  % analog beamforming
            end
            
            if ~obj.isLinear
                x = obj.elem{7}.step(x);  % Mixer thermal noise
                x = obj.elem{8}.step(x);  % Mixer AM/AM and AM/PM
            end
            
            % This stage is the AGC
            if obj.nbits ~= 0
                x = obj.elem{9}.step(x);
            end
            
            % The last stage of the RFFE is the ADC
            y = obj.elem{10}.step(x);
        end
    end
end