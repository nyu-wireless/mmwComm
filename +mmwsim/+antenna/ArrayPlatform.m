classdef ArrayPlatform < matlab.System
    % ArrayWithAxes.  Class containing an antenna array and axes.
    properties
        
        fc = 28e9;  % Carrier frequency
        
        % Element within each array.  Empty indicates to use an
        % isotropic element.  If non-empty, the object must be derived
        % from the matlab.System class with a step method
        %    dir = elem.step(fc,az,el)
        % that provides the directivity in dBi as a function of the
        % frequency and angles (az,el)
        elem = [];
        
        % Antenna array.  Empty indicates that there is only one element
        arr = [];
        
        % Steering vector object
        sv = []; 
        
        % Azimuth and elevation angle of the element peak directivity
        axesAz = 0;
        axesEl = 0;
        
        % Axes of the element local coordinate frame of reference
        axesLoc = eye(3);
        
        % Velocity vector in 3D in m/s
        vel = zeros(1,3);
        
        % Position in m
        pos = zeros(1,3);
    end
    
    methods
        function obj = ArrayPlatform(varargin)
            % Constructor
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
                
        function alignAxes(obj,az,el)
            % Aligns the axes to given az and el angles
            
            % Set the axesAz and axesEl to az and el
            obj.axesAz = az;
            obj.axesEl = el;
            
            % Creates axes aligned with az and el
            obj.axesLoc = azelaxes(az,el);
        end
        
        function dop = doppler(obj,az,el)
            % Computes the Doppler shift of a set of paths
            % The angles of the paths are given as (az,el) pairs
            % in the global frame of reference.
            
            % Finds unit vectors in the direction of each path
            npath = length(el);
            [u1,u2,u3] = sph2cart(deg2rad(az),deg2rad(el),ones(1,npath));
            u = [u1; u2; u3];
            
            % Compute the Doppler shift of each path via an inner product
            % of the path direction and velocity vector.
            vcos = obj.vel*u;
            vc = physconst('lightspeed');
            dop = vcos*obj.fc/vc;            
            
        end
        
        function releaseSV(obj)
            % Creates the steering vector object if it has not yet been
            % created.  Otherwise release it.  This is needed since the 
            % sv object requires that it takes the same number of 
            % inputs each time.
            if isempty(obj.sv)
                obj.sv = phased.SteeringVector('SensorArray',obj.arr);
            else
                obj.sv.release();
            end
            
        end
        
        function elemPosGlob = getElementPos(obj)
            % Gets the array elements in the global reference frame
           
            % Get the element position in the local reference frame
            elemPosLoc = obj.arr.getElementPosition();
            
            % Convert to the global reference frame
            elemPosGlob = local2globalcoord(elemPosLoc, 'rr', ...
                zeros(3,1), obj.axesLoc) + reshape(obj.pos,3,1);
        end
        
        function gain = getResponse(obj,az,el,w,relSV)
            % Computes the complex gain given angles and BF directions
            
            % Get default values
            if nargin < 4
                w = 1;
            end
            if nargin < 5
                relSV = false;
            end
            
            % Get SV and element gain
            [u, elemGain] = obj.step(az,el,relSV);
            
            % Compute gain
            elemGain = 10.^(0.05*elemGain);
            gain = (w.'*u) .* elemGain;
                        
            
        end
        
        function [gain,az,el] = getResponse2D(obj,az,el,w)
            % Gets the complex gain on a 2D angular grid
            % 
            % If w is nant x 1, this produces a gain matrix of 
            % size nel x naz representing the complex gain in 
            % each azimuth and elevation angle.  If w is nant x nw, 
            % gain is nw x nel x naz representing the complex gain in 
            % each azimuth and elevation angle and beamforming direction.
            
            % Set default angles to test
            if nargin < 2
                az = (-180:2:180)';
            end
            if nargin < 3
                el = (-90:2:90)';
            end
            if nargin < 4
                w = 1;
            end
            
            % Get grid of values
            [azMat, elMat] = meshgrid(az, el);
            azVal = azMat(:)';
            elVal = elMat(:)';
            
            % Get SV and element gain
            gain = obj.getResponse(azVal,elVal,w,true);
            
            % Reshape to a matrix
            nel = length(el);
            naz = length(az);
            nw = size(w,2);
            if (nw == 1)
                gain = reshape(gain, nel, naz);
            else 
                gain = reshape(gain, nw, nel, naz);
            end
        end
    end
    
    methods (Access = protected)
        
        function setupImpl(obj)
            % setup:  This is called before the first step.
            
            % Create the steering vector object using the array.
            obj.sv = phased.SteeringVector('SensorArray', obj.arr);
        end
        
        function releaseImpl(obj)
            % release:  Called to release the object
            obj.elem.release();
            obj.sv.release();
        end
        
       function [u, elemGain] = stepImpl(obj, az, el, relSV)
            % Gets steering vectors and element gains for a set of angles
            % The angles az and el should be row vectors along which
            % the outputs are to be computed.  
            % If the relSV == true, then the steering vector object is
            % released.  This is needed in case the dimensions of the past
            % call are the different from the past one
            
            % Release the SV
            if nargin < 4
                relSV = false;
            end
            if relSV
                obj.releaseSV();
            end
            
            % Convert the global angles (az, el) to local
            % angles (azLoc, elLoc).  Use the 
            % global2localcoord() method with the 'ss' option.
            uglobal = [az; el; ones(1,length(az))];
            ulocal = global2localcoord(uglobal, 'ss', zeros(3,1), obj.axesLoc);
            azLoc = ulocal(1,:);
            elLoc = ulocal(2,:);
            
            % TODO: Get the SV in the local coordinates
            %    u = obj.sv(...)
            u = obj.sv(obj.fc, [azLoc; elLoc]);  
                        
            % TODO:  Get the directivity gain of the element from the
            % local angles.
            %    elemGain = obj.elem.step(...) 
            elemGain = obj.elem.step(azLoc, elLoc);
        end

    end
    
    
end

