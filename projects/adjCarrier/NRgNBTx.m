classdef NRgNBTx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;	% carrier configuration
        pdschConfig;	% default PDSCH configuration
        waveformConfig;	% waveform configuration
        ofdmGridLayer;	% before pre-coding nsc x nsym x nlayers
        
        % Transmitted data
        bits;		% tx data bits
        pdschSym;	% tx data symbols
        dmrsSym;	% demodulation reference signals (DM-RS)
        ptrsSym;	% phase tracking reference signals (PT-RS)
        
        Nslot = 0;	% slot number
        fsamp;		% sample frequency of the DAC
        ncc;		% number of component carriers
        fcc;		% component carrier center frequency
        wtx;		% tx beamforming vector
        isInterference;
        Pnoise=1;
    end
    
    methods
        function obj = NRgNBTx(simParam, varargin)
            % Constructor
            
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the
            % objects need to be converted to older structures.
            obj.carrierConfig = mmwsim.nr.objtostruct(simParam.carrierConfig);
            obj.pdschConfig = mmwsim.nr.objtostruct(simParam.pdschConfig);
            obj.waveformConfig = mmwsim.nr.objtostruct(simParam.waveformConfig);
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % set the sample frequency.
            obj.fsamp = simParam.waveformConfig.SamplingRate;
            
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
        end
    end
    
    methods (Access = protected)
        
        function [x, z] = stepImpl(obj)
            % Set the slot number in the PDSCH config
            obj.pdschConfig.Nslot = obj.Nslot;
            
            % Get information for PDSCH and DM-RS allocations
            [pdschIndices, dmrsIndices, obj.dmrsSym, ...
                ptrsIndices, obj.ptrsSym, ~] = ...
                mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
                                   
            % Create the PDSCH grid before pre-coding
            obj.ofdmGridLayer = zeros(...
                obj.waveformConfig.NSubcarriers,...
                obj.waveformConfig.SymbolsPerSlot, ...
                obj.pdschConfig.NLayers);
            
            % Generate random bits
            bitsPerSym = mmwsim.nr.NRConst.bitsPerSym(...
                obj.pdschConfig.Modulation);
            nsym = length(pdschIndices);
            nbits = bitsPerSym * nsym;
            obj.bits = randi([0 1], nbits, 1);
            
             % Modulate the bits to symbols
            M = mmwsim.nr.NRConst.modOrder(obj.pdschConfig.Modulation);
            obj.pdschSym = qammod(obj.bits,M,'InputType','bit',...
                'UnitAveragePower',true);
                
            % Map symbols to OFDM grid
            obj.ofdmGridLayer(pdschIndices) = obj.pdschSym;
            obj.ofdmGridLayer(dmrsIndices) = obj.dmrsSym;
            obj.ofdmGridLayer(ptrsIndices) = obj.ptrsSym;
            
            % Perform the OFDM modulation
            xlayer = mmwsim.nr.hOFDMModulate(obj.carrierConfig, ...
                obj.ofdmGridLayer);
            
            % Upsample the Tx signal to the total bandwidth
            xup = resample(xlayer, obj.ncc, 1);
                
            % Use a NCO to shift the component carrier
            xnco = mmwsim.nr.hCarrierAggregationModulate(xup, ...
                obj.fsamp, obj.fcc(1));
            
            % Generate the interference z~CN(0,1)
            noiseGridLayer = sqrt(0.5)*randn(size(obj.ofdmGridLayer)) + 1j*randn(size(obj.ofdmGridLayer));
            zlayer = mmwsim.nr.hOFDMModulate(obj.carrierConfig, noiseGridLayer);
            zup = resample(zlayer, obj.ncc, 1);
            znco = mmwsim.nr.hCarrierAggregationModulate(zup, obj.fsamp, obj.fcc(2));
            
            % Remove mean
            xnco = xnco - mean(xnco);
            znco = znco - mean(znco);
            
            % Perform the TX beamforming.
            x = xnco*obj.wtx';
            z = znco;
            
            % Increment the slot number
            obj.Nslot = obj.Nslot + 1;
        end
    end
end