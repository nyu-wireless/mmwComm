classdef AqnmSim < matlab.System
    % aqnmSim: main symulation class
    
    properties
        % Simulation parameters
        noiseTemp = 290;	% noise temperature in K
        nrx = 16;			% number of RX antennas
        nx = 1e4;			% number samples per SNR point
        ktVar = 1.0;		% thermal noise per sample
        EkT;				% thermal noise energy
        xvar = 1;			% variance of the tx symbols
        snrInTest;			% input SNR
        
        % RFFE
        rffe;				% rf front-end object of the receiver
        isLinear;			% include the rffe non-linearities
        lnaGain;			% lna gain in dB
        lnaNF;				% lna noise figure in dB
        lnaAmpLut;			% lna fund tone power curve
        lnaPower;			% lna power in mW
        
        mixGain;			% mixer gain in dB
        mixNF;				% mixer noise figure in dB
        mixAmpLut;			% mixer fund tone power curve
        mixPLO;				% mixer input power from the local oscillator
        mixPower;			% mixer power in mW
        
        % ADC
        fsamp;				% sample frequency of the ADC
        adcFOM = 65e-12;	% FOM of the ADC 65 fJ/conv
        nbits = 4;			% ADC bits per dimension
        
        % Transmit symbol type:  'iidGaussian' or 'iidPhase'
        txSymType = 'iidGaussian';
        
        % Channel type: 'iidGaussian', 'iidPhase', 'iidAoA' or 'ones'
        chanType = 'iidPhase';
    end
    
    methods
        function obj = AqnmSim(varargin)
            % Constructor
            
            % Set parameters from constructor arguments.
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            obj.EkT = physconst('Boltzman')*obj.noiseTemp;
            
            % Create the RX RFFE
            obj.rffe = mmwsim.rffe.RFFERx(...
                'nrx', obj.nrx, ...
                'lnaNF', obj.lnaNF, ...
                'lnaGain', obj.lnaGain, ...
                'lnaPower', obj.lnaPower, ...
                'lnaAmpLut', obj.lnaAmpLut, ...
                'mixNF', obj.mixNF, ...
                'mixPLO', obj.mixPLO, ...
                'mixGain', obj.mixGain, ...
                'mixPower', obj.mixPower, ...
                'mixAmpLut', obj.mixAmpLut, ...
                'fsamp', obj.fsamp, ...
                'nbits', obj.nbits, ...
                'isLinear', obj.isLinear);
        end
        
        function [y,w,x] = genData(obj)
            % Generate the random channel
            if strcmp(obj.chanType,'iidPhase')
                phase = 2*pi*rand(obj.nx,obj.nrx);
                w = exp(1i*phase);
            elseif strcmp(obj.chanType, 'iidGaussian')
                w = (randn(obj.nx,obj.nrx) + 1i*randn(obj.nx,obj.nrx))/sqrt(2);
            elseif strcmp(obj.chanType, 'randAoA')
                dsep = 0.5;
                theta = unifrnd(-pi/2,pi/2,obj.nx,1);
                phase = 2*pi*cos(theta)*(0:obj.nrx-1)*dsep;
                w = exp(1i*phase);
            elseif strcmp(obj.chanType, 'ones')
                w = ones(obj.nx,obj.nrx);
            else
                error('Unknown channel type');
            end
            
            % Generate random symbols without scaling
            if strcmp(obj.txSymType, 'iidGaussian')
                x = (randn(obj.nx,1) + 1i*randn(obj.nx,1))*sqrt(1/2);
            elseif strcmp(obj.txSymType, 'iidPhase')
                x = exp(1i*2*pi*rand(obj.nx,1));
            else
                error('Unknown TX symbol type');
            end
            
            % Generate RX symbols with no noise
            y = x.*w;
            
            % Rescale so that it is Es/kT = 1
            scale = sqrt(obj.EkT/mean(abs(y).^2, 'all'));
            y = y * scale;
        end
        
        function NF = nf(obj)
            % Calculate the effective noise figure
            NF = obj.rffe.nf();
        end
        
        function P = power(obj)
            P = obj.rffe.power();
        end
    end
    
    methods (Access = protected)
        function [snrOut] = stepImpl(obj)
            [y, w, x] = obj.genData();
            
            % Find the number of rf elements
            nsnr = length(obj.snrInTest);
            snrOut = zeros(nsnr, 1);
            
            for isnr = 1:nsnr
                % Get the SNR and scale the input signal
                ynoisy = 10^(0.05*obj.snrInTest(isnr))*y;
                
                % Run through RFFE stages
                ynoisy = obj.rffe.step(ynoisy);
                
                % Beamform
                xhat = sum(ynoisy.*conj(w),2) ./ sum(abs(w).^2,2);
                
                % Use a linear receiver
                %
                % xhat = a*x + d,  d ~ CN(0, E|xhat-x|^2)
                a = mean(conj(xhat).*x)/mean(abs(x).^2);
                dvar = mean(abs(xhat - a*x).^2);
                
                % Measure the output SNR
                snrOut(isnr) = 10*log10(abs(a).^2*obj.xvar/dvar);
                
                % Display progress
                if false
                    fprintf(1, 'snrIn = %12.4e, snrOut = %12.4e\n', ...
                        obj.snrInTest(isnr), snrOut(isnr));
                end
            end
        end
    end
end

