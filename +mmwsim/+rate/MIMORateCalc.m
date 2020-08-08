classdef MIMORateCalc < matlab.System
    % MIMORateCalc Class for computing the narrowband MIMO rate
    %
    % The class computes the rate for a narrowband model of the form
    %    y = Hx + w
    % The model assumes perfect CSIR and CSIT (H is known at both TX and
    % RX).  In this case, it performs and SVD and then computes the rate on
    % each virtual stream.
    %
    properties
        nstreamsMax = 4;  % maximum number of streams
        
        % Includes rate loss from Shannon limit.
        includeRateLoss = true;
        
        % If includeRateLoss == True, the losses are described by
        % the following parameters
        specEffMax = 4.8; % maximum spectral efficiency per stream
        snrLoss = 3;      % SNR loss from Shannon limit
        bwLoss = 0.2;     % bandwidth overhead
                
        % Parameters for last capacity calculation
        U, V;       % Linear TX and RX matrices
        s;          % singular values along H
        rateStr;    % Rate as a function of the number of streams
        
    end
    
    methods
        function obj = MIMORateCalc(varargin)
            % Constructor:
            % The syntax allows you to call the constructor with syntax of
            % the form:
            %
            %     chan = MIMOMPChan('Prop1', Val1, 'Prop2', val2, ...);
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
        
        function rate = computeRate(obj,H,snr)
            % Computes the rate of the channel. 
            %
            % The snr=Es/N0 in dB where
            % Es = E\|x\|^2 is the TX symbol energy across all TX antennas
            % and N0 = E|w(k)|^2 is the RX noise energy per antenna

            % Take the SVD
            [obj.U, S, obj.V] = svd(H, 0);
            obj.s = diag(S);
            
            % Compute the SNR in linear domain
            if obj.includeRateLoss
                snrLin = 10^(0.1*(snr-obj.snrLoss));
            else
                snrLin = 10^(0.1*snr);
            end
                            
            % Compute rate as a function of the number of streams
            if obj.nstreamsMax > 0
                r = min(length(obj.s), obj.nstreamsMax);            
            else
                r = length(obj.s);
            end
            obj.rateStr = zeros(r,1);            
            for i = 1:r
                % Capacity assuming power is allocated uniformly
                % across all streams
                c = log2(1 + obj.s(1:i).^2*snrLin/i);
                if obj.includeRateLoss
                    c = min(c, obj.specEffMax);
                    obj.rateStr(i) = (1-obj.bwLoss)*sum(c);
                else
                    obj.rateStr(i) = sum(c);
                end                                                                    
            end
            
            rate = max(obj.rateStr);
        end
    end
end

