classdef NRUERx < matlab.System
    % 5G NR UE receiver class
    
    properties
        % Configuration
        carrierConfig;  % Carrier configuration
        pdschConfig;    % PDSCH configuation
        waveformConfig; % Waveform configuration
        
        % RFFE
        nstages = 1;
        rffe;     % rf front-end object
        phy;      % physical-layer object
        ps;       % phase shifter model
        lna;      % LNA model
        mixer;    % Mixer model
        nbits;    % ADC num of bits
        isLinear; % include the rffe non-linearities
        isFD;     % fully-digital or analog beamforming
        nrx;
    end
    
    methods (Access = public)
        
        % Constructor
        function obj = NRUERx(varargin)
            % Constructor
            
            % Set the parameters from the constructor arguments
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            % Create the RX RFFE
            obj.rffe = mmwsim.rffe.RFFERx(...
                'nrx', obj.nrx, ...
                'ps', obj.ps, ...
                'lna', obj.lna, ...
                'nstages', obj.nstages, ...
                'mixer', obj.mixer, ...
                'SampleRate', obj.waveformConfig.SampleRate, ...
                'nbits', obj.nbits, ...
                'isLinear', obj.isLinear, ...
                'isFD', obj.isFD);
            
            % Create the RX PHY
            obj.phy = mmwsim.phy.PHYRx(...
                'ncc', obj.waveformConfig.ncc, ...
                'isFixPoint', false, ...
                'enableCA', true, ...
                'NRB', obj.carrierConfig.NSizeGrid, ...
                'SCS', obj.carrierConfig.SubcarrierSpacing, ...
                'fsamp', obj.waveformConfig.SampleRate, ....
                'nbadc', obj.nbits, ...
                'fcc', -obj.waveformConfig.fcc);
            
        end
        function pdschEq = stepImpl(obj, rxWaveform, wrx)
            % Demodulate one slot of data
            
            % Find the DMRS and PTRS symbols and indices.
            dmrsSym = nrPDSCHDMRS(obj.carrierConfig, obj.pdschConfig);
            dmrsInd = nrPDSCHDMRSIndices(obj.carrierConfig, obj.pdschConfig);
            ptrsSym = nrPDSCHPTRS(obj.carrierConfig, obj.pdschConfig);
            ptrsInd = nrPDSCHPTRSIndices(obj.carrierConfig, obj.pdschConfig);
            
            [pdschInd, pdschInfo] = nrPDSCHIndices(obj.carrierConfig, obj.pdschConfig);
            
            % Pass the data from the RFFE layer
            if obj.isFD
                rxWaveformRFFE = obj.rffe.step(rxWaveform) * wrx;
            else
                rxWaveformRFFE = obj.rffe.step(rxWaveform, wrx);                
            end
            
            % Pass the data from the PHY layer
            rxWaveformPHY = obj.phy.step(rxWaveformRFFE);
            sig = rxWaveformPHY(:,1);
            
            % Demodulate the RX signal
            rxGrid = nrOFDMDemodulate(obj.carrierConfig, sig);
            
            % Estimate channel and noise variance
            [chanEstGrid, noiseEst] = nrChannelEstimate( ...
                obj.carrierConfig, rxGrid, dmrsInd, dmrsSym, ...
                'CDMLengths', obj.pdschConfig.DMRS.CDMLengths);
            
            % Create a temporary grid
            tempGrid = nrResourceGrid(obj.carrierConfig, obj.pdschConfig.NumLayers);
            
            % MMSE Equalization of PDSCH
            [pdschRx, pdschHest] = nrExtractResources(pdschInd, rxGrid, chanEstGrid);
            pdschEq = nrEqualizeMMSE(pdschRx, pdschHest, noiseEst);
            tempGrid(pdschInd) = pdschEq;
            
            % MMSE Equalization of PTRS
            [ptrsRx, ptrsHest] = nrExtractResources(ptrsInd, rxGrid, chanEstGrid);
            ptrsEq = nrEqualizeMMSE(ptrsRx, ptrsHest, noiseEst);
            tempGrid(ptrsInd) = ptrsEq;
            
            % Estimate the Common Phase Error (CPE)
            cpe = nrChannelEstimate(obj.carrierConfig, tempGrid, ptrsInd, ptrsSym);
            cpe = angle(sum(cpe, [1 3 4]));
            
            % Compensate CPE
            symLoc = pdschInfo.PTRSSymbolSet(1)+1:pdschInfo.PTRSSymbolSet(end)+1;
            tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1j*cpe(symLoc));
            
            % Update the equalized symbols
            pdschEq = tempGrid(pdschInd);
        end
    end
end