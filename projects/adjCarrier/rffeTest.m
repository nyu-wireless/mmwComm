%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
fc = 140e9; % carrier frequency in Hz
noiseTemp = 290; % noise temperature in K
EkT = physconst('Boltzman')*noiseTemp; % noise energy
nslot = 20; % number of slots
ndsn = 2;
nrx = 1; % num of receive antennas
wtx = 1; % TX beamforming vector
wrx = 1; % RX beamforming vector
isLinear = false;
isFixPoint = false;
enableCA = true;
isFilter = false;

% Parameters describing the PDSCH
simParam = PDSCHSimParam('fc', fc);
fsamp = simParam.waveformConfig.SamplingRate;
ncc = simParam.ncc;

% Since, many 5G Toolbox routines do not take classes, the objects need to
% be converted to older structures.
carrierConfig = mmwsim.nr.objtostruct(simParam.carrierConfig);
pdschConfig = mmwsim.nr.objtostruct(simParam.pdschConfig);
waveformConfig = mmwsim.nr.objtostruct(simParam.waveformConfig);

% Load the nonlinear models
load('rffe140GHz.mat');

%% Main Loop
adcTest = [4, 5, 6]; % ADC resolution
snrSigTest = linspace(-10, 50, 7)'; % Signal power in dB
snrIntTest = linspace(-10, 50, 7)'; % Interference in dB

nsig = length(snrSigTest);
nint = length(snrIntTest);
snrOut = zeros(nsig, nint, ndsn, nslot);

% Create a TX object using the NRgNBTx class
tx = NRgNBTx(simParam, ...
    'wtx', wtx, ...
    'ncc', ncc);

% Create the receiver objects using the NRUERx class
rx = cell(ndsn, 1);

for idsn = 1:ndsn
    if idsn == 1
        ilna = 4;
        imix = 5;
        iplo = 2;
        iadc = 1;
    else
        ilna = 3;
        imix = 5;
        iplo = 7;
        iadc = 2;
    end
    
    rx{idsn} = NRUERx(...
        'nrx', nrx, ...
        'ncc', ncc, ...
        'wrx', wrx, ...
        'lnaNF', lnaNF(ilna), ...
        'lnaGain', lnaGain(ilna), ...
        'lnaPower', lnaPower(ilna), ...
        'lnaAmpLut', lnaAmpLut(:,:,ilna), ...
        'mixNF', mixNF(iplo), ...
        'mixPLO', mixPLO(iplo), ...
        'mixGain', mixGain(imix), ...
        'mixPower', mixPower(imix), ...
        'mixAmpLut', reshape(mixAmpLut(:,:,iplo,imix),31,3), ...
        'fsamp', fsamp, ...
        'isLinear', isLinear, ...
        'isFixPoint', isFixPoint, ...
        'enableCA', enableCA, ...
        'nbits', adcTest(iadc), ...
        'waveformConfig', waveformConfig, ...
        'carrierConfig', carrierConfig, ...
        'pdschConfig', pdschConfig, ...
        'isFilter', isFilter);
end

for islot = 1:nslot
    fprintf(1, '%2d: \t', islot);
    tic;
    % Generate the TX sequence and the interference
    [x, z] = tx.step();
    
    % Rescale so that it is Es/kT = 1
    x = x * sqrt((fsamp*EkT)/mean(abs(x).^2, 'all'));
    z = z * sqrt((fsamp*EkT)/mean(abs(z).^2, 'all'));
    
    for idsn = 1:ndsn
        rx{idsn}.pdschSymTx = tx.pdschSym;
        
        for isig = 1:nsig
            % Scale the signal power
            xnoisy = x .* 10^(0.05*snrSigTest(isig));
            for iint = 1:nint
                % Scale the noise power
                znoisy = z .* 10^(0.05*snrIntTest(iint));
                
                % The output SNR at that signal and interference level
                snrOut(isig, iint, idsn, islot) = rx{idsn}.step(xnoisy+znoisy);
            end
        end
    end
    toc;
end

% Take the average over multiple slots
snrOut = mean(snrOut, 4);

% Save the noise figure for each receiver design
nf = [rx{1}.nf(), rx{2}.nf()];

%% Clean up
clear adcTest carrierConfig enableCA fc fsamp iadc idsn iint ilna imix iplo;
clear isFilter isFixPoint isig isLinear islot lnaAmpLut lnaFOM lnaGain;
clear lnaIIP3 lnaNF EkT noiseTemp;
clear lnaPower mixAmpLut mixGain mixIIP3 mixNF mixPLO mixPower ncc ndsn;
clear nint nrx nsig nslot pdschConfig rnd rx simParam tx waveformConfig;
clear wrx wtx x xnoisy z znoisy
