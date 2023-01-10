classdef NRgNBTx < matlab.System
    % 5G NR gNB transmitter class
    
    % Public properties
    properties
        % Configuration
        carrierConfig;  % Carrier configuration
        pdschConfig;  % PDSCH configuation
        waveformConfig;  % Waveform configuration
        pdschSym;
        fcc;
    end
    
    methods (Access = public)
        % Constructor
        function obj = NRgNBTx(varargin)
            % Set the parameters from the constructor arguments
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function [txWaveform] = stepImpl(obj)
            % 
            % wtx: precoding matrix
            % Find the DMRS and PTRS symbols and indices.
            dmrsSym = nrPDSCHDMRS(obj.carrierConfig, obj.pdschConfig);
            dmrsInd = nrPDSCHDMRSIndices(obj.carrierConfig, obj.pdschConfig);
            ptrsSym = nrPDSCHPTRS(obj.carrierConfig, obj.pdschConfig);
            ptrsInd = nrPDSCHPTRSIndices(obj.carrierConfig, obj.pdschConfig);
            
            % Create an OFDM grid
            txGrid = nrResourceGrid(obj.carrierConfig, ...
                obj.pdschConfig.NumLayers);
            
            % Get the indices where the PDSCH is allocated
            [pdschInd, pdschInfo] = nrPDSCHIndices(obj.carrierConfig, ...
                obj.pdschConfig);
            
            % Generate random bits
            txBits = randi([0 1], pdschInfo.G, 1);
            
            % Modulate the bits to symbols
            obj.pdschSym = nrPDSCH(obj.carrierConfig, obj.pdschConfig, txBits);
            
            % Map the modulated symbols to the OFDM grid
            txGrid(pdschInd) = obj.pdschSym;
            txGrid(dmrsInd) = dmrsSym;
            txGrid(ptrsInd) = ptrsSym;
            
            % Modulate
            txWaveform = nrOFDMModulate(obj.carrierConfig, txGrid);
            
            % Increase the slot number
            obj.carrierConfig.NSlot = obj.carrierConfig.NSlot + 1;
        end
    end
end