classdef DAC < matlab.System
    % DAC class with scaling.
    properties	
        % number of bits, 0 indicates no quantization
        nbits = 6;
        isComplex = true;    % complex input                 
        phaseDither = false; % Enable phase dithering

        % DAC scaling
        aq = 1.0;			% quantizer step
        nscal = 1000;       % number of samples used for calibration

        % Parameters for linear model:
        %    Q(x) = linGain*x + N(0,quantVar),   x~N(0,inputVar)
        inputVar = 1;       % input variance        
        linGain = 1;    
        quantVar = 1;
        mseOpt = 0;         % optimal MSE in dB
        qvarSim;

        fsamp = 491.52e6;   % sample frequency
    end

    % Quantizer methods.  All methods are static
    methods 
        function obj = DAC(varargin)
        % Constructor
        % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end      
        end

        % Quantizer with saturation
        % Quantizes input x, with number of bits, nb and scale a.
        function y = qsat(obj, x0)
            % Scale the input        
            x = x0 / sqrt(obj.inputVar);
            if (obj.nbits == 0)
                % No quantization
                y = x;
            elseif ~obj.isComplex
                % Perform quantization for real signals
                y = (max( min(round(obj.aq*x-0.5),2^(obj.nbits-1)-1),-2^(obj.nbits-1))+0.5)/obj.aq;
            else
                % If phase dither is enabled, we rotate the input 
                % by random phases
                if obj.phaseDither
                    phase = 2*pi*rand(size(x));
                    x = x.*exp(1i*phase);
                end

                % Perform for IQ
                xr = real(x);
                xi = imag(x);
                yr = (max( min(round(obj.aq*xr-0.5),2^(obj.nbits-1)-1),-2^(obj.nbits-1))+0.5)/obj.aq;
                yi = (max( min(round(obj.aq*xi-0.5),2^(obj.nbits-1)-1),-2^(obj.nbits-1))+0.5)/obj.aq;                
                y = yr + 1i*yi;                

                % Unwrap the phase dithering
                if obj.phaseDither
                    y = y.*exp(-1i*phase);
                end                               
            end
            % Rescale with input variance
            y = y * sqrt(obj.inputVar);                                                                                  
            obj.qvarSim = mean(abs(y-x0).^2,2);
        end

        % Sets scaling values from a different DAC
        function copyScale(obj, adc1)
            obj.aq = adc1.aq;
            obj.inputVar = adc1.inputVar;
            obj.quantVar = adc1.quantVar;
            obj.linGain = adc1.linGain;
            obj.mseOpt = adc1.mseOpt;
        end

        % Finds the optimal quantizer scale level.
        % The optimal value is found by searching over values aq to
        % minimize
        %   mse = E(x - qsat(x,nb,aq))^2
        % when x = N(0,inputVar).  
        function optScale(obj)
            % Handle case with infinite resolution
            if (obj.nbits == 0)
                obj.aq = 1;
                obj.quantVar = 0;
                obj.linGain = 1;
                obj.mseOpt = 0;
                return
            end

            % Generate random data points to calibrate
            if obj.isComplex
                x = randn(1, obj.nscal) + 1i*randn(1,obj.nscal);
                x = x * sqrt(obj.inputVar/2);
            else
                x = randn(1, obj.nscal) * sqrt(obj.inputVar);
            end

            % Measure MSE on possible test quantizer levels
            aqtest = linspace(0.1, 2, 500)'*2^(obj.nbits-1);
            naq  = length(aqtest);
            mse = zeros(naq,1);

            for i = 1:naq
                obj.aq = aqtest(i);
                mse(i) = mean( abs(x-obj.qsat(x)).^2 );
            end

            % Select scaling with minimal distortion
            [~, im] = min(mse);  
            obj.aq = aqtest(im);

            % Find parameters for a linear AQN model
            %    qsat(x) = linGain*x + N(0, quantVar)                     
            q = obj.qsat(x);
            obj.linGain  = real(q*x')/real(x*x');
            obj.quantVar = mean(abs(q - obj.linGain*x).^2 );

            % Compute the relative MMSE:
            %    mseOpt = var(x|q)/var(x)
            obj.mseOpt = 10*log10( obj.quantVar/...
            (obj.inputVar*obj.linGain^2 + obj.quantVar) );
        end
    end

    methods (Access = protected)
        function y = stepImpl(obj, x)
            % Step function:  Performs the quantization
            y = obj.qsat(x);
        end
    end
end