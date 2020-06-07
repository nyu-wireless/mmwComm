classdef InterpPatternAntenna < matlab.System
    % Wrapper class for an antenna to perform smooth interpolation of the
    % directivity.  This class is necessary since the antenna toolbox
    % generally performs nearest neighbor interpolation which may not
    % be smooth
    properties
        
        ant;  % Base antenna.  Must support a pattern method
        fc;   % Frequency        
        dirInterp;  % Gridded interplant
    end
    
    methods
        function obj = InterpPatternAntenna(ant, fc, varargin)
            % Constructor
            obj.ant = ant;
            obj.fc = fc;
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
    end
        
    methods (Access = protected)
        
        function setupImpl(obj)
            % setup:  This is called before the first step.            
            % We will use this point to create the interpolator.
            
            % Get the pattern from ant.pattern
            [dirPat,azPat,elPat] = obj.ant.pattern(...
                obj.fc,'Type', 'Directivity');                     
            
            % Create the gridded interpolant object.
            obj.dirInterp = griddedInterpolant({elPat,azPat},dirPat);         
            
        end
        function dir = stepImpl(obj, az, el)
            % Computes the directivity along az and el angles in the local
            % reference frame         
            
            % Run the interplation object to compute the directivity
            % in the local angles
            dir = obj.dirInterp(el,az);            
        end
    end
end

