classdef Quant < matlab.System
    % Quantizer methods.  All methods are static
    methods (Static)
        
        % Quantizer with saturation
        % quantizes input x, with number of bits, nb and scale a.
        function y = qsat(x, nb, aq)
            if (nargin < 3)
                aq = 1;
            end
            if (nb == 0)
                % No quantization
                y = x;
            elseif isreal(x)
                % Perform quantization for real signals
                y = (max( min(round(aq*x-0.5),2^(nb-1)-1),-2^(nb-1))+0.5)/aq;
            else
                % Perform for IQ
                xr = real(x);
                xi = imag(x);
                yr = (max( min(round(aq*xr-0.5),2^(nb-1)-1),-2^(nb-1))+0.5)/aq;
                yi = (max( min(round(aq*xi-0.5),2^(nb-1)-1),-2^(nb-1))+0.5)/aq;
                y = yr + 1i*yi;                
            end
        end
        
        % Finds the optimal quantizer scale level
        % The optimal value is found by searching over values aq to
        % minimize
        %   mse = E(x - qsat(x,nb,aq))^2
        % when x = N(0,1).  
        function [aqopt, msemin] = optScale(nb, ns)
            % Generate random Gaussian vector
            if (nargin < 2)
                ns = 1e4;
            end
            x = randn(1, ns);
            
            % Measure MSE on possible test quantizer levels
            aqtest = linspace(0.1, 2, 500)'*2^(nb-1);
            naq  = length(aqtest);
            mse = zeros(naq,1);
            for i = 1:naq
                aq = aqtest(i);
                mse(i) = mean( abs(x-mmwsim.rffe.Quant.qsat(x, nb, aq)).^2 );
            end
            
            % Select scaling with minimal distortion
            [msemin, im] = min(mse);            
            aqopt = aqtest(im);
        end
    end
end