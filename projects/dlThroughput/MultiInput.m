classdef MultiInput <  matlab.System
	% MultiInput:  Replicates a MATLAB system object with multiple inputs
	properties
		nin = 1;  % number of inputs

		% If clone==true, then sys is a cell array of systems,
		% with sys{i} being the system for input i
		% If clone==false, then the system is applied to all the inputs
		% Use this only when the systems can share the same state,
		% or the system is stateless.
		clone = true;
		sys;
	end

	methods
		
		function obj = MultiInput(sys0,nin,clone)
			% Constructor
			%
			% Clones the system, sys0, for all the inputs.
			
			if nargin >= 2
				obj.nin = nin;
			end
			
			if nargin >= 3
				obj.clone = clone;
			end
			
			if obj.clone
				obj.sys = cell(obj.nin,1);
				for i = 1:obj.nin
					obj.sys{i} = sys0.clone();
				end
			else
				obj.sys = obj.sys0;
			end
		end
	end


	methods (Access = protected)

		function y = stepImpl(obj, x)
			% Step function:  Performs the system on each input

			% Check input dimensions
			if size(x,2) ~= obj.nin
				error('Expected %d inputs received %d',...
				obj.nin, size(x,2));
			end

			% Apply system to each input
			nt = size(x,1);
			y = zeros(nt,obj.nin);
			
			for i = 1:obj.nin
				if obj.clone
					sysi = obj.sys{i};
				else
					sysi = obj.sys;
				end
				y(:,i) = sysi.step(x(:,i));
			end
		end
	end
end

