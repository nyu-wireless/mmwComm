classdef ADC < matlab.System
    % ADC class with scaling.
    
    properties
        % number of bits, 0 indicates no quantization
        nbits = 6;
        
        % Output type:
        % "int":  signed int of the form,
        %     q = 2*k+1, k=-M/2 to M/2-1, M=2^nbits
        % "float":  floating point output of the form:
        %     q = stepSize*(k+0.5)
        % If isComplex, then this representation is used in I and Q
        outputType = "int";
        
        
        isComplex = true;    % complex input
        stepSize = 1.0;          % step size
        dither = false;      % enable dithering
        
        % ADC scaling parameters
        nscal = 10000;       % number of samples used for calibration
        
        % Parameters for linear model:
        %    Q(x) = linGain*x + N(0,quantVar),   x~N(0,inputVar)
        inputVar = 1;       % input variance
        linGain = 1;
        quantVar = 1;
        quantVar1 = 1;
        mseOpt = 0;         % Optimal MSE in dB
    end
    
    % Quantizer methods.  All methods are static
    methods
        
        function obj = ADC(varargin)
            % Constructor
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
        
        % Performs the quantization providing integer and floating
        % point values
        function [qfloat,qint] = qsat(obj, x0)
            if (obj.nbits == 0)
                % No quantization
                qint = x0;
                qfloat = x0;
                return
            end
            
            
            % Scale the input
            x = x0 / obj.stepSize;
            
            % Dithering
            if obj.dither
                if obj.isComplex
                    %d = (rand(size(x))-0.5) + 1i*(2*rand(size(x))-0.5);
                    d = (rand(size(x))) + 1i*(2*rand(size(x)));
                else
                    d = rand(size(x))-0.5;
                end
                x = x + d;
            end
            
            M2 = 2^(obj.nbits-1);
            if obj.isComplex
                % Perform quantization for complex signals
                xr = floor(real(x));
                xi = floor(imag(x));
                qr = max(min(xr,M2-1), -M2)+0.5;
                qi = max(min(xi,M2-1), -M2)+0.5;
                qint = qr + 1i*qi;
            else
                % Perform quantization for complex signals
                x = floor(real(x));
                qint = max(min(x,M2-1), -M2)+0.5;
            end
            
            % Scale output to integer or float
            qint = 2*qint;
            qfloat = 0.5*obj.stepSize*qint;
            
            % Remove dithering
            if obj.dither
                qfloat = qfloat - d*obj.stepSize;
            end
            
        end
        
        % Sets scaling values from a different ADC
        function copyScale(obj, adc1)
            obj.stepSize = adc1.stepSize;
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
                obj.stepSize = 1;
                obj.quantVar = 0;
                obj.linGain = 1;
                obj.mseOpt = 0;
                return
            end
            
            % Generate random data points to calibrate
            if obj.isComplex
                x = randn(1, obj.nscal) + 1i*randn(1,obj.nscal);
                xstd = sqrt(obj.inputVar/2);
                x = x * xstd;
                
            else
                xstd = sqrt(obj.inputVar);
                x = randn(1, obj.nscal) * xstd;
            end
            
            % Measure MSE on possible test quantizer levels
            stepTest = linspace(0, 4, 500)'*xstd/2^(obj.nbits-1);
            ns  = length(stepTest);
            mse = zeros(ns,1);
            for i = 1:ns
                obj.stepSize = stepTest(i);
                mse(i) = mean( abs(x-obj.qsat(x)).^2 );
            end
            
            % Select scaling with minimal distortion
            [~, im] = min(mse);
            obj.stepSize = stepTest(im);
            
            % Find parameters for a linear AQN model
            %    qsat(x) = linGain*x + N(0, quantVar)
            [qf, qi] = obj.qsat(x);
            if strcmp(obj.outputType, 'int')
                q = qi;
            else
                q = qf;
            end
            obj.linGain  = real(q*x')/real(x*x');
            obj.quantVar = mean(abs(q - obj.linGain*x).^2 );
            obj.quantVar1 = mean(abs(q - x).^2 );
            
            % Compute the relative MMSE:
            %    mseOpt = var(x|q)/var(x)
            obj.mseOpt = 10*log10( obj.quantVar/...
                (obj.inputVar*obj.linGain^2 + obj.quantVar) );
            
        end
    end
    
    methods (Access = protected)
        
        function q = stepImpl(obj, x)
            % Step function:  Performs the quantization
            [qf, qi] = obj.qsat(x);
            if strcmp(obj.outputType, 'int')
                q = qi;
            else
                q = qf;
            end
        end
    end
end