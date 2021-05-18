%% 5G NR downlink PDSCH - Non Linear RFFE
% This project evaluates a practical mobile receiver at millimeter-wave
% (mmWave) and terahertz (THz) frequencies. We consider a downlink system
% with a single NR basestation (gNB), and a single UE device.
% The simulation contains realistic models of RF front-end (RFFE)
% components for the receiver.
%
% The gNB can either use the entire wideband bandwidth by performing
% carrier aggregation, or transmit a single component carrier. The gNB
% generates for each component carrier a physical downlink shared channel
% (PDSCH) that includes the information data and some physical layer
% signals. The demodulation reference signals (DM-RS) and the phase
% tracking reference signals (PT-RS). The UE uses DM-RS and PT-RS to
% perform practical channel estimation. For each slot we generate a random
% angle-of-departure (AoD) and angle-of-arrival (AoA) in both azimuth and
% elevation angles. For each RFFE configuration we create a UE device that
% will process the received data at several input power levels.

%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
% We will use the following parameters
fc = 140e9;			% carrier frequency in Hz
noiseTemp = 290;	% noise temperature in K
EkT = physconst('Boltzman')*noiseTemp; % noise energy

nantUE = [4,4];		% array size at the UE (mobile device)
nantgNB = [8,8];	% array size at the gNB (base station)
ueVel = [5 0 0];	% UE velocity vector in m/s
nrx = nantUE(1)*nantUE(2);	% num of receive antenna

isFixPoint = false;	% use low resolution PHY-layer processing
isLOS = true;		% use a single line-of-sight path
enableCA = false;	% enable carrier aggregation
isHPC = false;		% run the simulation at the NYU HPC

% Input SNR Es/N0 relative to thermal noise
snrInTest = linspace(-10, 70, 21)';

% ADC resolution (4-bit, 5-bit, 6-bit). For inf-bit use 0.
adcTest = [4, 5, 6];

% Define a uniform random function in [a,b]
rnd = @ (a,b) (a + (b-a)*rand(1,1));

% Load the RFFE models
if fc == 140e9
    load('rffe140GHz.mat');
    isLinear = false;	% 'false' to include the distortion from the rffe
else
    % if there are no models for this carrier frequency use ideal
    % components
    isLinear = true;
end

% Parameters describing the PDSCH
simParam = PDSCHSimParam('fc', fc, 'enableCA', enableCA);
fsamp = simParam.waveformConfig.SamplingRate;
ncc = simParam.ncc;

% Load a NR channel
dlySpread = 50e-9;  % delay spread in seconds
chan = nrCDLChannel('DelayProfile','CDL-A',...
    'DelaySpread',dlySpread, 'CarrierFrequency', fc, ...
    'NormalizePathGains', true);
chaninfo = info(chan);

% Get parameters from simulation parameters. Since, many 5G Toolbox
% routines do not take classes, the objects need to be converted to older
% structures.
carrierConfig = mmwsim.nr.objtostruct(simParam.carrierConfig);
pdschConfig = mmwsim.nr.objtostruct(simParam.pdschConfig);
waveformConfig = mmwsim.nr.objtostruct(simParam.waveformConfig);

% Compute received input power in dBm
Pin = 10*log10(EkT*fsamp) + 30 + snrInTest;

if isHPC
    arrayId = getenv('SLURM_JOB_ID');
    if isempty(aID)
        arrayId = randi([0,100]);
    end
    rng(str2double(arrayId),'twister');
    nit = 20;	% number of slots
else
    rng('shuffle');
    arrayId = 0;
    nit = 5;	% number of slots
end

%% Create the antenna arrays for the gNB and the UE

% Constants
vp = physconst('lightspeed');	% speed of light
lambda = vp/fc;					% wavelength

% Create a patch element
len = 0.49*lambda;
groundPlaneLen = lambda;
elem = patchMicrostrip( ...
    'Length', len, 'Width', 1.5*len, ...
    'GroundPlaneLength', groundPlaneLen, ...
    'GroundPlaneWidth', groundPlaneLen, ...
    'Height', 0.01*lambda, ...
    'FeedOffset', [0.25*len 0]);

% Tilt the element so that the maximum energy is in the x-axis
elem.Tilt = 90;
elem.TiltAxis = [0 1 0];

% Create wrapper class
elemInterp = mmwsim.antenna.InterpPatternAntenna(elem, fc);

% The arrays are separated at lambda/2.
dsep = 0.5*lambda;
arrUE0  = phased.URA(nantUE,dsep,'ArrayNormal','x');
arrgNB0 = phased.URA(nantgNB,dsep,'ArrayNormal','x');

% Create a wrapper class around the arrays to handle orientation modeling.
% The ArrayPlatform class includes an array along with axes to store the
% orientation of the array relative to some global coordinate system.
arrgNB = mmwsim.antenna.ArrayPlatform('arr', arrgNB0, ...
    'elem', elemInterp, 'fc', fc);
arrUE  = mmwsim.antenna.ArrayPlatform('arr', arrUE0, ...
    'elem', elemInterp, 'fc', fc);
arrUE.set('vel', ueVel);

if ~isLOS
    % Get the gains and other path parameters
    gain = chaninfo.AveragePathGains;
    aoaAz  = chaninfo.AnglesAoA;
    aoaEl = 90-chaninfo.AnglesZoA;
    aodAz  = chaninfo.AnglesAoD;
    aodEl = 90-chaninfo.AnglesZoD;
    dly = chaninfo.PathDelays;
    
    % Find the paths with maximum gain
    [gainmax, im] = max(gain);
end
%% Main simulation loop
nsnr = length(snrInTest);			% num of snr tests
nadc = length(adcTest);				% num of adc
nlna = size(lnaAmpLut, 3);			% num of lna
nmix = size(mixAmpLut, 4);			% num of mixer
nplo = length(mixPLO);				% num of input LO power to the mixer
ndrivers = length(0:log2(nrx));		% num of LO driver configurations
nsim = nlna * nmix * nplo * nadc;	% num of simulations

% Intialize the output SNR vector
snrOut = zeros(nsnr, nit);

iadc = 1;
ilna = 1;
imix = 1;
iplo = 1;

% Loop over the slots
for it = 1:nit
    fprintf('\nSlot: %d\n',it);
    tic;
    
    if isLOS
        gain = 0;
        aoaAz = rnd(-180, 180);
        aoaEl = rnd(-90, 90);
        aodAz = rnd(-180, 180);
        aodEl = rnd(-90, 90);
        dly = 0;
        im = 1;
    end
    
    % Find the spatial signatures and element gains of each path based on
    % their angles of arrival and departure.
    [utx, elemGainTx] = arrgNB.step(aodAz, aodEl);
    [urx, elemGainRx] = arrUE.step(aoaAz, aoaEl);
    
    % Compute the TX beamforming direction at the gNB and RX beamforming
    % direction at the UE.  To keep the beamforming simple, we will align
    % the directions to the strongest path.  Thus, the BF directions should
    % be complex conjugate of the steering vectors.  They should also be
    % normalized.
    wtx = conj(utx(:,im)); % TX direction at the gNB
    wtx = wtx / norm(wtx);
    wrx = conj(urx(:,im)); % RX direction at the UE
    wrx = wrx / norm(wrx);
    
    % Create a TX object using the NRgNBTx class
    tx = NRgNBTx(simParam, 'wtx', wtx, 'ncc', ncc);
    
    % Create the channel
    chan = mmwsim.chan.MIMOMPChan('aoaAz', aoaAz, 'aoaEl', aoaEl, ...
        'aodAz', aodAz, 'aodEl', aodEl, 'gain', gain, 'rxArr', arrUE, ...
        'txArr', arrgNB, 'dly', dly, 'fsamp', fsamp);
    
    % Create TX and RX
    x = tx.step();
    y = chan.step(x);
    
    % Rescale so that it is Es/(fs*kT) = 1
    scale = sqrt((fsamp*EkT)/mean(abs(y).^2, 'all'));
    y = y * scale;
    
    % Create all possible receiver configurations
    rx = NRUERx(...
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
        'snrInTest', snrInTest, ...
        'pdschSymTx', tx.pdschSym, ...
        'waveformConfig', waveformConfig, ...
        'carrierConfig', carrierConfig, ...
        'pdschConfig', pdschConfig);
    
    % If it is possible run parrallel simulations to calculate the output
    % SNR of a linear receiver
    snrOut(:, it) = rx.step(y);	% Output SNR [dB]
    
    % print the processing time for each slot
    toc;
end

% Average over all iterations.
snrOut = mean(snrOut, 2);
%% Plot the Output SNR

figure(1); clf; 
plot(Pin, snrOut, '-', 'linewidth', 1.5);
box on;
axis tight;
xlabel('Receive power per antenna [dBm]', 'interpreter', 'latex', ...
	'fontsize', 13);
ylabel('Output SNR $\;(\gamma_\mathrm{out})\;$ [dB]', ...
	'interpreter', 'latex', ...
	'fontsize', 13);
grid on;