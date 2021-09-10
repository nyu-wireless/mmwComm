classdef IIRFilter < matlab.System
    %IIRFILTER
    
    properties
        fp = 190.08; % Passband Frequency
        fst = 245.76; % Stopband Frequency
        ap = 1; % Passband Ripple (dB)
        ast = 80; % Stopband Attenuation (dB)
        fs = 491.52; % Sampling Frequency
        hd;
    end
    
    methods
        function obj = IIRFilter(varargin)
            % Constructor
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
    end
    methods (Access = protected)
        
        function setupImpl(obj)
            
            h = fdesign.lowpass('fp,fst,ap,ast', ...
                obj.fp, obj.fst, obj.ap, obj.ast, obj.fs);
            
            obj.hd = design(h, 'ellip');
            
        end
        
        function y = stepImpl(obj, x)
            y = filter(obj.hd, x);
        end
    end
end

