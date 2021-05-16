classdef PHYRx < matlab.System
    % PhyRx. The physical layer processing after the A/D
    
    properties
        NRB;			% num of resource blocks
        SCS;			% subcarrier spacing
        isFixPoint;		% use low resolution PHY-layer processing
        
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
        enableCA;
        ncc;		% number of component carriers
        fcc;		% component carrier center frequency
        
        % ADC
        fsamp;
        nbadc;
    end
    
    methods
        function obj = PHYRx(varargin)
            % Constructor
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            % Find the center frequency for each component carrier.
            if obj.ncc == 4
                obj.fcc = [-150e6, -50e6, 50e6, 150e6];
            elseif obj.ncc == 8
                obj.fcc = [-700e6, -500e6, -300e6, -100e6, ...
                    100e6, 300e6, 500e6, 700e6];
            elseif obj.ncc == 2
                obj.fcc = [-100e6, 100e6];
            else
                obj.fcc = 0;
            end
            obj.fcc = -obj.fcc;
            
            % Create a fixed-point FIR low pass filter for low-resolution
            % baseband processing.
            if obj.isFixPoint
                % Find the effective signal bandwidth:
                % ResourceBlocks * 12 * SubCarrierSpacing
                fsig = obj.NRB * 12 * obj.SCS * 1e3;
                
                obj.Fp = fsig/obj.fsamp;	% pass-band frequency
                obj.Fst = 1/obj.ncc;		% stop-band frequency
                
                % Design a Fixed-Point Filter
                spec = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
                    obj.Fp, obj.Fst, obj.Ap, obj.Ast);
                f = design(spec, 'minphase', false, 'SystemObject', true);
                
                coef = coeffs(f);
                
                bq = fi(coef.Numerator, true, obj.nbcoeff, ...
                    'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
                L = bq.FractionLength;
                bsc = coef.Numerator*2^L;
                
                obj.hd = dfilt.dffir(bsc);
                obj.hd.Arithmetic = 'fixed';
                obj.hd.CoeffWordLength = obj.nbcoeff;
                
                % Integer real input from ADC with nbits resolution
                obj.hd.InputWordLength = obj.nbadc;
                obj.hd.InputFracLength = 0;
                
                obj.bfilt = obj.hd.Numerator;
                obj.afilt = 1;
            end
        end
        
        function y = stepImpl(obj, x)
            if obj.enableCA
                % reshape the input signal to allow seperate processing for
                % each component carrier.
                y = zeros(size(x,1)/obj.ncc, size(x,2), obj.ncc);
                
                for icc = 1:obj.ncc
                    % Use a NCO to bring the component carrier of interest
                    % to the middle of the spectrum. This, allows us to use
                    % the same optimized low-pass filter for all the
                    % component carriers.
                    xnco = mmwsim.nr.hCarrierAggregationModulate(...
                        x, obj.fsamp, obj.fcc(icc));
                    
                    if obj.isFixPoint
                        % filter and then downsample the input signal with
                        % a filter that has low-resolution coefficients
                        xfilt = filter(obj.bfilt, obj.afilt, xnco);
                        y(:,:,icc) = downsample(xfilt, obj.ncc);
                    else
                        % use an ideal low-pass filter to downsample the
                        % input signal.
                        y(:,:,icc) = sqrt(obj.ncc)*resample(xnco, 1, obj.ncc); %
                    end
                end
            else
                y = x;
            end
        end
    end
end