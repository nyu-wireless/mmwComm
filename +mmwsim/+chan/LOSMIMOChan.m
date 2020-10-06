classdef LOSMIMOChan < matlab.System
    % LOSMIMOChan:  High-rank narrowband MIMO channel
    properties 
        fsamp = 122e6;  % Sample rate in Hz
        fc = 28e9;      % Carrier frequency in Hz
        
        % TX and RX antenna arrays
        txArr, rxArr;
        
        % Fixed gain.  If not specified, then the gain is computed
        % from Friis' law and the element gains
        gainFix = [];  % Fixed gain
        
		fracDly;
		gain = 0;
		
        % Latest channel parameters
        chanMatrix;  % Narrowband channel matrix
        dly;         % Delay in seconds
        gainTx, gainRx;  % elementa gains in dBi
        aodAz, aodEl, aoaAz, aoaEl;  % Angles of arrival and departure
        pathLoss;  % free space omni path loss in dB
        distCen;   % distance between array centers in m
        utx, urx;  % Spatial signatures (these are not needed)                       
    end
    
    methods 
        function obj = LOSMIMOChan(varargin)
            % Constructor:  
            % The syntax allows you to call the constructor with syntax of
            % the form:
            %
            %     chan = MIMOMPChan('Prop1', Val1, 'Prop2', val2, ...);
            if nargin >= 1
                obj.set(varargin{:});
            end
            
        end
        
        function computeChanMatrix(obj, fcResp)
            % Computes the narrowband channel matrix and delay
            % The matrix is stored in chanMatrix
            
            % Default parameters
            if nargin < 2
                fcResp = obj.fc;
            end
            
            % Get the element positions in the global coordinate system
            txpos = obj.txArr.getElementPos();
            rxpos = obj.rxArr.getElementPos();
            
            % Compute the distances between the elements
            ntx = size(txpos,2);
            nrx = size(rxpos,2);
            d = reshape(rxpos,3,nrx,1) - reshape(txpos,3,1,ntx);
            dist = sqrt(sum(d.^2,1));
            dist = reshape(dist,nrx,ntx);

            % Get the phase change along each path
            vc = physconst('lightspeed');
            lambda = vc/fcResp;
            phase = 2*pi*dist/lambda;
            
            if isempty(obj.gainFix)
                % Gain is not fixed.  Then we compute delay and gain
                % from Friis' law
                       
                % Get the separation vector from the centers of the arrays
                % We will use this single value for the free space path loss
                % and element responses
                txcen = mean(txpos,2);
                rxcen = mean(rxpos,2);
                dcen = rxcen - txcen;

                % Convert to spherical coordinates
                [obj.aodAz, obj.aodEl, obj.distCen] = ...
                    cart2sph(dcen(1), dcen(2), dcen(3));
                [obj.aoaAz, obj.aoaEl, ~] = ...
                    cart2sph(-dcen(1), -dcen(2), -dcen(3));  
                obj.aodAz = rad2deg(obj.aodAz);
                obj.aodEl = rad2deg(obj.aodEl);
                obj.aoaAz = rad2deg(obj.aoaAz);
                obj.aoaEl = rad2deg(obj.aoaEl);

                % Get the TX steering vectors and element gains along the 
                % angles of departure using the        
                [obj.utx, obj.gainTx] = obj.txArr.step(obj.aodAz, obj.aodEl);            
                [obj.urx, obj.gainRx] = obj.rxArr.step(obj.aoaAz, obj.aoaEl);

                % Compute free space path loss from the center direction
                obj.pathLoss = fspl(obj.distCen, lambda);

                % Create the narrowband channel matrix
                gainLin = 10^(0.05*(obj.gainTx+obj.gainRx-obj.pathLoss));
                obj.chanMatrix = gainLin*exp(1i*phase); 

                % Compute the absolute delay
                obj.dly = obj.distCen / vc;
            else
                
                % Gain is fixed.
                gainLin = 10.^(0.05*obj.gainFix);
                obj.dly = 0;
            end
			
            % Compute the channel matrix
            obj.chanMatrix = gainLin*exp(1i*phase);
                        
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