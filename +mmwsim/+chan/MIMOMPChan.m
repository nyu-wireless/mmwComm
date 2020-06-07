classdef MIMOMPChan < matlab.System
    % MIMOMPChan:  MIMO multi-path fading channel    
    properties 
        fsamp;   % Sample rate in Hz
        
        % TX and RX antenna arrays
        txArr, rxArr;
                
        % Path properties
        aoaAz, aoaEl; % Angles of arrival in degrees
        aodAz, aodEl; % Angles of departure in degrees
        gain;  % path gains in dB
        dly;   % delays in seconds
        dop;   % doppler shift of each path in Hz
        
        % Fractional delay object
        fracDly;
        
        % Initial set of phases for the next step call
        phaseInit;
        
        % Previous gains and steering vectors
        utx, urx;
        gainTx, gainRx;
        
                                
    end
    
    methods 
        function obj = MIMOMPChan(varargin)
            % Constructor:  
            % The syntax allows you to call the constructor with syntax of
            % the form:
            %
            %     chan = MIMOMPChan('Prop1', Val1, 'Prop2', val2, ...);
            if nargin >= 1
                obj.set(varargin{:});
            end
            
        end
        
    end
    methods (Access = protected)
        function setupImpl(obj)
              % setup:  This is called before the first step.
              
              % Create a dsp.VariableFractionalDelay object 
              obj.fracDly = dsp.VariableFractionalDelay(...
                'InterpolationMethod', 'Farrow','FilterLength',8,...
                'FarrowSmallDelayAction','Use off-centered kernel',...
                'MaximumDelay', 1024);                           
        end
        
        function resetImpl(obj)
            % reset:  Called on the first step after reset or release.
            
            % Reset the fracDly object
            obj.fracDly.reset();
            
            % Initialize phases, phaseInit, to a row vector of 
            % dimension equal to the number of paths with uniform values 
            % from 0 to 2pi
            npath = length(obj.gain);
            obj.phaseInit = 2*pi*rand(1,npath);
        end
        
        function releaseImpl(obj)
            % release:  Called after the release method
            
            % Release the fracDly object
            obj.fracDly.release();
        end
        
        function y = stepImpl(obj, x)
            % step:  Run samples through the channel
            % The input, x, should be nsamp x nanttx, 
            
                        
            % Compute the delay in samples
            dlySamp = obj.dly*obj.fsamp;                        
            
            % Get the TX steering vectors and element gains along the 
			% angles of departure using the        
            [obj.utx, obj.gainTx] = obj.txArr.step(obj.aodAz, obj.aodEl);            
            [obj.urx, obj.gainRx] = obj.rxArr.step(obj.aoaAz, obj.aoaEl);
            
            % TODO:  Compute the total gain along each path in linear scale
            % The gain is the path gain + element gains at the TX and RX
            %    gainLin = ...
            gainLin = 10.^(0.05*(obj.gain + obj.gainTx + obj.gainRx));
            
                        
            % Initialize variables    
            nsamp  = size(x,1);
            nantrx = size(obj.urx,1);
            npath = length(obj.dly);
            y = zeros(nsamp,nantrx);
            
            % Get the Doppler shift of each path from the TX and RX using 
			% txArr.doppler() and rxArr.doppler() methods
            obj.dop = obj.rxArr.doppler(obj.aoaAz, obj.aoaEl) + ...
                      obj.txArr.doppler(obj.aodAz, obj.aodEl);
            
            % Using the Doppler shifts, compute the phase rotations 
            % on each path.  Specifically, if nsamp = length(x), create a
            % (nsamp+1) x npath matrix 
            %     phase(i,k) = phase rotation on sample i and path k
            nsamp = size(x,1);
            phase = 2*pi*(0:nsamp)'*obj.dop/obj.fsamp + obj.phaseInit;
            
            % Save the final phase, phase(nsamp+1,:) as phaseInit for the 
			% next step.
            obj.phaseInit = phase(nsamp+1,:);
            
            % Loop over the paths
            for ipath = 1:npath

                % Compute the transmitted signal, x, along the 
                % TX spatial signature for path ipath.
                z = x*obj.utx(:,ipath);

                % Delay the path by the dlySamp(ipath) using the fractional
				% delay object
                z = obj.fracDly(z,dlySamp(ipath));
                
                % Multiply by the gain 
                z = gainLin(ipath)*z;
                
                % Rotate by the phase 
                z = z .* exp(1i*phase(1:nsamp,ipath));
                
                % Multiply by the RX spatial signature and add to y
                y = y +  z*obj.urx(:,ipath).';
            end
          
        end
    end
end