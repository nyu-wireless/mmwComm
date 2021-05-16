classdef Codebook < matlab.mixin.SetGetExactNames
    % Codebook:  Class for generating codebooks for multi-array
    % antenna systems
    properties
        % Test directions
        ntest = 5000;  % number of test vectors
        d = 3;         % dimension of the space
        
        % Nuumber of iterations
        nit = 20;
        
        % Test directions wrt universal coordinates
        X0;             % Cartesian
        azTest, elTest; % azimuth and elevation angles
        
        % Array placement parameters
        % Cell array of arrays
        arrSet = {};
        narr;  % number of arrays
        
        % Codebook design parameters
        nant;       % number of antennas in each array.
        ncode = []; % number of codewords per array
        Z;  % Z(i,j,k) = array resp from dir i antenna j array k
        W;  % W(j,ell,k) = codebook ell
        
        % Codebook computed values
        rhoMax;  % rhoMax(i) = highest correlation for dir i
        kmax;    % kmax(i) = best array index
        jmax;    % jmax(i,k) = best codeword for dir i, array k
    end
    
    methods
        function obj = Codebook(varargin)
            % Constructor
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
        
        function genTestDir(obj)
            % Generates test directions to evaluate beamforming
            
            % Generate test vectors
            % The test vectors are used to compute the codebook
            X = randn(obj.ntest,obj.d);
            Xnorm = sqrt(sum(abs(X).^2,2));
            obj.X0 = X ./ repmat(Xnorm,[1,obj.d]);
            
            % Compute angles
            [obj.azTest, obj.elTest, ~] = ...
                cart2sph(obj.X0(:,1), obj.X0(:,2), obj.X0(:,3));
            obj.azTest = rad2deg(obj.azTest);
            obj.elTest = rad2deg(obj.elTest);
            
        end
        
        
        
        function genCodebook(obj)
            % Generates the codebook
            
            % Get the number of antennas
            obj.narr = length(obj.arrSet);
            for i = 1:obj.narr
                nanti = obj.arrSet{1}.arr.getNumElements();
                if i == 1
                    obj.nant = nanti;
                elseif (obj.nant ~= nanti)
                    e = MException('Codebook',...
                        'Currently all arrays must have same num elements');
                    throw(e);
                end
            end
            
            % If not specified, set number of codewords = number of
            % antennas
            if isempty(obj.ncode)
                obj.ncode = obj.nant;
            end
            
            % Generate vectors to test
            obj.genTestDir();
            
            % Get the array response from each array
            %   Z(itest,ielem,iarr) = phase and element gain
            %      from direction itest, element ielem and array iarr
            obj.Z = zeros(obj.ntest, obj.nant, obj.narr);
            for i = 1:obj.narr
                
                [Ui,elemGaini] = obj.arrSet{i}.step(...
                    obj.azTest', obj.elTest',true);
                obj.Z(:,:,i) = (Ui.') .* 10.^(0.05*elemGaini');
            end
            
            % Compute the initial codebook
            obj.W = zeros(obj.nant, obj.ncode, obj.narr);
            for i = 1:obj.narr
                % Generate a random codebook
                Wi = randn(obj.nant, obj.ncode) + 1i*randn(obj.nant, obj.ncode);
                Wnorm = sqrt(sum(abs(Wi).^2,1));
                obj.W(:,:,i) = Wi ./ repmat(Wnorm,[obj.nant,1]);
            end
            
            % Run k-means
            
            for it = (1:obj.nit)
                
                % Compute correlations per array
                %  rho(i,j,k) = corr for dir i, code j on array k
                %  jmax(i,k) = best codeword for direction i on array k
                %  rhoMax0(i,k) = best corr for direction i on array k
                rho = zeros(obj.ntest,obj.ncode,obj.narr);
                obj.jmax = zeros(obj.ntest,obj.narr);
                rhoMax0 = zeros(obj.ntest,obj.narr);
                for k = 1:obj.narr
                    rho(:,:,k) = abs(obj.Z(:,:,k)*obj.W(:,:,k)).^2;
                    [rhoMax0(:,k), obj.jmax(:,k)] = max(rho(:,:,k),[],2);
                end
                
                % Compute best across arrays
                [obj.rhoMax, obj.kmax] = max(rhoMax0, [], 2);
                
                % Update the codewords for each array
                for k = 1:obj.narr
                    % Directions assigned to array k
                    Ik = find(obj.kmax == k);
                    nk = length(Ik);
                    
                    for j = 1:obj.ncode
                        % Directions assigned to array k and code j
                        Ijk = find((obj.kmax == k) & (obj.jmax(:,k)==j));
                        
                        if isempty(Ijk)
                            % If there are no directions assigned to the
                            % codeword, select a random direction assigned
                            % to the array
                            i = randi(nk,1,1);
                            wj = obj.Z(i,:,k)';
                            
                        else
                            % Otherwise, assign the codeword to the maximum direction
                            % for the set that it is assigned to
                            Zj = obj.Z(Ijk,:,k);
                            
                            % First singular vector of Zj gives direction of max
                            % correlation
                            [~, ~, wj] = svds(Zj,1);
                            
                        end
                        
                        % Re-assign the normalized codeword
                        wj = wj/sqrt(sum(abs(wj).^2));
                        obj.W(:,j,k) = wj;
                    end
                    
                end
                
            end
            
        end
        
        function [bfGainMax, bfGainArr, indCode, indArr] = getBfGain(obj,az,el)
            % Gets the array and codewords for a set of AoAs
            
            % Optimal gain and codeword index per array
            ndir = length(az);
            indCode = zeros(ndir,obj.narr);
            bfGainArr = zeros(ndir,obj.narr);
            
            for i = 1:obj.narr
                
                % Get the array response = sv + element gain
                [Ui,elemGaini] = obj.arrSet{i}.step(az, el, true);
                Zi = (Ui.') .* 10.^(0.05*elemGaini');
                
                % Get BF gain on each codeword
                rho = abs(Zi(:,:)*obj.W(:,:,i)).^2;
                
                % Find the best codeword within the array
                [bfGainLin, indCode(:,i)] = max(rho,[],2);
                bfGainArr(:, i)  = 10*log10(bfGainLin);
            end
            
            % Get optimal over all arrays
            [bfGainMax, indArr] = max(bfGainArr,[],2);
        end
    end
end