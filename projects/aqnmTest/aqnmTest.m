%% RFFE with Non-Linear Components
% This project evaluates a single-input multi-output (SIMO) system at
% millimeter-wave (mmWave) and terahertz (THz) frequencies. The simulation
% contains realistic models of RF front-end (RFFE) components for the
% receiver. For every RFFE configuration we measure the end-to-end
% performance and the receiver power consumption. We develop a generic
% model for characterizing the performance using two key parameters.

%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
fc = 140e9;			% carrier frequency in Hz
nrx = 16;			% number of RX antennas
isLinear = false;	% include the RFFE non-linear distortion

% Input SNR Es/N0 relative to thermal noise
snrInTest = linspace(-10, 50, 31)';

% ADC resolution (4-bit, 5-bit, 6-bit). For inf-bit use 0.
adcTest = [4, 5, 6];

% Transmit symbol type:  'iidGaussian' or 'iidPhase'
txSymType = 'iidGaussian';

% Channel type: 'iidGaussian', 'iidPhase', 'iidAoA' or 'ones'
chanType = 'iidPhase';

% Load the RFFE models
if fc == 140e9
    load('rffe140GHz.mat');
    
    % We consider a receiver at 140 GHz for 6G communications. The
    % bandwidth is calculated based on 3GPP NR waveforms with carrier
    % aggregation.
    fsamp = 1.96608e+09;
else
    % if there are no models for this carrier frequency use ideal
    % components
    isLinear = true;
end

% Compute received input power Pin in dBm
noiseTemp = 290;						% noise temperature in K
EkT = physconst('Boltzman')*noiseTemp;	% noise energy
Pin = 10*log10(EkT*fsamp) + 30 + snrInTest;


isHPC = false;	% 'true' to run the simulation at the NYU HPC

if isHPC
    arrayId = getenv('SLURM_JOB_ID');
    if isempty(aID)
        arrayId = randi([0,100]);
    end
    rng(str2double(arrayId),'twister');
    nit = 20;	% number of iterations
else
    rng('shuffle');
    arrayId = 0;
    nit = 5;	% number of iterations
end

%% Main simulation loop
nsnr = length(snrInTest);			% num of snr tests
nadc = length(adcTest);				% num of adc
nlna = size(lnaAmpLut, 3);			% num of lna
nmix = size(mixAmpLut, 4);			% num of mixer
nplo = length(mixPLO);				% num of input LO power to the mixer
ndrivers = length(0:log2(nrx));		% num of LO driver configurations
nsim = nlna * nmix * nplo * nadc;	% num of simulations

% Intialize vectors
snrOut = zeros(nsnr, nsim, nit);
sim = cell(nsim, 1);

for it = 1:nit
    fprintf('\nIteration: %d\n',it);
    tic;
    % Create all possible combinations.
    isim = 1;
    for iadc = 1:nadc
        for ilna = 1:nlna
            for imix = 1:nmix
                for iplo = 1:nplo
                    sim{isim} = AqnmSim(...
                        'nrx', nrx, ...
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
                        'nbits', adcTest(iadc), ...
                        'txSymType', txSymType, ...
                        'chanType', chanType, ...
                        'snrInTest', snrInTest);
                    isim = isim + 1;
                end
            end
        end
    end
    
    % If it is possible run parrallel simulations using `parfor`
    parfor isim = 1:nsim
        snrOut(:, isim, it) = sim{isim}.step(); % Output SNR [dB]
    end
    
    % print the processing time for each slot
    toc;
end

% Average over all iterations.
snrOut = mean(snrOut, 3);

% Calculate the saturation SNR
snrSat = snrOut(end, :)';

%% Find the Effective Noise Figure and Power Consumption

rffePower = zeros(nsim, ndrivers);
rffeNF = zeros(nsim, 1);

for isim = 1:nsim
    rffePower(isim, :) = sim{isim}.power();	% RFFE power consumption [mW]
    rffeNF(isim) = sim{isim}.nf();			% Effective noise figure [dBm]
end

% Find the minimum power for each parameter setting. The `idriver` will
% denote the number of LO drivers.
[rffePower, idriver] = min(rffePower, [], 2);

%% Fit a model
% We use two parameters to characterize the performance of each
% configuration: (a) effective noise figure that is dominant in low input
% power; (b) saturation SNR that is dominant in high SNR. Using these
% values we can fit a model for the output SNR as follows,

% For non-linear systems we empirically show that the output SNR can be
% calculated by the following formula
nom = nrx * 10.^(0.1*snrInTest);
denom = reshape(10.^(0.1*rffeNF), [], nsim) + ...
    nrx * reshape(10.^(0.1*snrInTest), nsnr, []) .* ...
    reshape(10.^(-0.1*snrSat), [], nsim);
rffeModel = 10*log10(nom./denom);

%% Save the data

% clear some variables
clear sim it hpc isim isnr ilna imix iplo iadc nom denom noiseTemp EkT;

% save the simulation results
save(sprintf('aqnmData%d.mat', arrayId));