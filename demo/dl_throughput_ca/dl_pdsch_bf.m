%% Demo: 5G NR Downlink Simulation with Beamforming 
%
% Beamforming is an essential component of the millimeter wave (mmWave)
% communication systems. This demo will demonstrate simple beamforming
% and channel modeling for a downlink transmissions in the 5G New Radio
% standard.
%% Clear the command window, the workspace and close all figures
clc;		% clear the command window
clear;		% clear all variables and functions from memory
close all;	% close all figures

%% Packages
%
% The repository contains a folder |+mmwsim| containing some routines 
% for simulating 5G NR systems.

% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
%
% We will use the following parameters
fc = 28e9;          % carrier frequency in Hz
nantUE = [4,4];     % array size at the UE (mobile device)
nantgNB = [8,8];    % array size at the gNB (base station)
snrPerAntenna = -5; % target SNR per sample per antenna  
ueVel = [5 0 0];    % UE velocity vector in m/s
subcarrierSpacing = 120;  % sub-carrier spacing in kHz
numCC = 4;	% number of component carriers

% Creates simulation parameters for this demo
simParam = PDSCHSimParam('fc', fc);

%% Load the 3GPP NR channel model
%
% We will load the same channel as in the previous lab.  
dlySpread = 50e-9;  % delay spread in seconds
chan = nrCDLChannel('DelayProfile','CDL-A',...
    'DelaySpread',dlySpread, 'CarrierFrequency', fc, ...
    'NormalizePathGains', true);
chaninfo = info(chan);

% Get the gains and other path parameters
gain = chaninfo.AveragePathGains;
aoaAz  = chaninfo.AnglesAoA;
aoaEl = 90-chaninfo.AnglesZoA;
aodAz  = chaninfo.AnglesAoD;
aodEl = 90-chaninfo.AnglesZoD;
dly = chaninfo.PathDelays;

im = 2;
gain = gain(im);
aoaAz = aoaAz(im);
aoaEl = aoaEl(im);
aodAz = aodAz(im);
aodEl = aodEl(im);
dly = dly(im);

%% Create the antenna arrays for the gNB and the UE

% Constants
vp = physconst('lightspeed');  % speed of light
lambda = vp/fc;   % wavelength

% Create a patch element
len = 0.49*lambda;
groundPlaneLen = lambda;
elem = patchMicrostrip(...
    'Length', len, 'Width', 1.5*len, ...
    'GroundPlaneLength', groundPlaneLen, ...
    'GroundPlaneWidth', groundPlaneLen, ...
    'Height', 0.01*lambda, ...
    'FeedOffset', [0.25*len 0]);

% Tilt the element so that the maximum energy is in the x-axis
elem.Tilt = 90;
elem.TiltAxis = [0 1 0];

% A problem with the Phased Array Toolbox is that the element patterns
% are not smoothly interpolated.  Also, we want to support antennas
% that have analytic functions for their pattern.  To support this, 
% we will use a wrapper class, InterpPatternAntenna.  This class
% derives from the system.Matlab super-class.  Its step method provides
% the directivity as a function of the angles.  We can then replace
% this class with any other class that provides a formula for the
% directivity.
elemInterp = mmwsim.antenna.InterpPatternAntenna(elem, fc);

% The arrays are separated at lambda/2.
dsep = 0.5*lambda;
arrUE0  = phased.URA(nantUE,dsep,'ArrayNormal','x');
arrgNB0 = phased.URA(nantgNB,dsep,'ArrayNormal','x');

% Create a wrapper class around the arrays to handle orientation modeling.
% The ArrayPlatform class includes an array along with axes to store the 
% orientation of the array relative to some global coordinate system.
arrgNB = mmwsim.antenna.ArrayPlatform('arr', arrgNB0, 'elem', elemInterp, 'fc', fc);
arrUE  = mmwsim.antenna.ArrayPlatform('arr', arrUE0, 'elem', elemInterp, 'fc', fc);
arrUE.set('vel', ueVel);

% Rotate the UE and gNB antennas to the path with the maximum gain. Find 
% the index of the path with the maximum gain.
[gainmax, im] = max(gain);

% Call the arrUE.alignAxes() and arrgNB.alignAxes() to  align to the 
% corresponding angles of arrival and departure.
arrUE.alignAxes(aoaAz(im), aoaEl(im));
arrgNB.alignAxes(aodAz(im), aodEl(im));

%% Compute the gains along the paths
% To see the potential gain from beamforming, we will compute the 
% array factor and element gain along each path

% Find the spatial signatures and element gains of each path based on their
% angles of arrival and departure.
[utx, elemGainTx] = arrgNB.step(aodAz, aodEl);
[urx, elemGainRx] = arrUE.step(aoaAz, aoaEl);

% Compute the TX beamforming direction at the gNB and RX BF direction at 
% the UE.  To keep the beamforming simple, we will align the directions to 
% the strongest path.  Thus, the BF directions should be complex conjugate
% of the steering vectors.  They should also be normalized.
wtx = conj(utx(:,im)); % TX direction at the gNB
wtx = wtx / norm(wtx);
wrx = conj(urx(:,im)); % RX direction at the UE
wrx = wrx / norm(wrx);

% Compute the array factors at the gNB and UE from the BF vectors and 
% spatial signatures, utx and urx.
%    AFgNB(i) = array factor gain on path i in dBi at the gNB
%    AFUE(i) = array factor gain in path i dBi at the UE
AFgNB = 20*log10(abs(wtx.'*utx));
AFUE = 20*log10(abs(wrx.'*urx));

% Compute the gain on each path adding the array factors and elemement 
% gains.
gainDir = gain + AFgNB + AFUE + elemGainTx + elemGainRx;

% Plot both the original gain and gainDir, the gain with directivity.
% stem(dly/1e-9, [gain; gainDir]', 'BaseValue', -40);
% grid on;
% xlabel('Delay (ns)');
% ylabel('Gain (dB)');
% legend('Omni', 'With directivity');
%% Generate a 5G TX signal 
%
% We will now test the array processing by transmitting a 5G downlink
% signal.  Specifically, we will transmit random QPSK symbols on the
% locations of the PDSCH channel, the channel in the 5G NR standard
% for data.  Most of the class is implemented and extensively uses
% commands from the 5G Toolbox.

% Create a TX object using the NRgNBTx object for each CC. 
tx = NRgNBTx(simParam);

% Set the BF vector of the TX
tx.set('txBF', wtx);

% Generate one slot of symbols
x = tx.step();

specPlot = mmwsim.nr.hCarrierAggregationPlotSpectrum(x(:,1), tx.fsamp,...
    'Power Spectrum of Carrier Aggregation Waveform',...
    {'Carrier aggregated signal'});

%% Create the MIMO multi-path channel
% We will now simulate the channel in time-domain.  The lab supplies
% code, MIMOMPChan, which is a MIMO version of the SISOMPChan
% you created in the previous lab.

% Create the MIMOMPChan object with all the AoAs, AoDs, gains arrays, 
% delays and sampling rate.  
%    chan = ...
chan = mmwsim.chan.MIMOMPChan('aoaAz', aoaAz, 'aoaEl', aoaEl, 'aodAz', aodAz, ...
    'aodEl', aodEl, 'gain', gain, 'rxArr', arrUE, 'txArr', arrgNB, ...
    'dly', dly, 'fsamp', tx.fsamp);

y = chan.step(x);

%% Add noise
%
% In multi-antenna receivers, the SNR is typically quoted as the SNR per
% antenna.  Specifically, suppose that
%    ynoisy = y + w.

% The SNR per antenna is E|y(t,j)|^2/E|w(t,j)|^2.
Erx = mean(abs(y).^2, 'all');
Enoise = 10.^(-0.1*snrPerAntenna)*Erx;
ynoisy = y + sqrt(Enoise/2)*(randn(size(y)) +1i*randn(size(y)));
%% Create a UE receiver
%
% We will now demodulate the noisy symbols.   The lab supplies a simple 
% class, NRUErx, to perform this function.  Most of the class is 
% implemented and extensively uses commands from the 5G Toolbox.  You just 
% have to do a small modification to support BF.  

% Create a RX object with the correct rxBF vector
rx = NRUERx(simParam, 'rxBF', wrx);

% Run the rx.step() method with ynoisy to receive the signals
rx.step(ynoisy);

% The equalized symbols are now stored in rx.pdschSym. 
subplot(2,2,1);
plot(rx.pdschSymEq(:, 1), '.');
title('Component Carrier = 1');
xlabel('real');
ylabel('imag');
xlim([-2,2]);
ylim([-2,2]);
grid();

subplot(2,2,2);
plot(rx.pdschSymEq(:, 2), '.');
title('Component Carrier = 2');
xlabel('real');
ylabel('imag');
xlim([-2,2]);
ylim([-2,2]);
grid();

subplot(2,2,3);
plot(rx.pdschSymEq(:, 3), '.');
title('Component Carrier = 3');
xlabel('real');
ylabel('imag');
xlim([-2,2]);
ylim([-2,2]);
grid();

subplot(2,2,4);
plot(rx.pdschSymEq(:, 4), '.');
title('Component Carrier = 4');
xlabel('real');
ylabel('imag');
xlim([-2,2]);
ylim([-2,2]);
grid();


% Measure the SNR
% When you plot the final equalized symbols you will see that there
% is a lot of noise.  Also there is a phase rotation which comes
% from the Doppler shift that was not corrected.  In reality, you would
% have some carrier frequency offset to remove this.  We will also discuss
% this later.  For now, we measure the post-equalized SNR.
%
% One way to measure the post-equalized SNR is:
%    snrEq = 10*log10( E|r|^2 / E| r - h*x |^2 )
% where r is the recived raw symbols (in this case rx.pdschSymRaw)
% h is the channel estimate (rx.pdschChanEst) and 
% x is the transmitted symbols (tx.pdschSym).

% Compute and print the post-equalized SNR in dB
for ncc = 1:4
	Eerr = mean(abs(rx.pdschSymRaw(:, ncc) - rx.pdschChanEst(:, ncc).*tx.pdschSym(:, ncc)).^2);
	Erx = mean(abs(rx.pdschSymRaw(:, ncc)).^2);
	snrEq = 10*log10(Erx / Eerr);
	fprintf(1,'SNR post equalization = %7.2f\n', snrEq);
end

nantrx = prod(nantUE);
bwGain = 10*log10(tx.waveformConfig.Nfft/tx.waveformConfig.NSubcarriers);
snrTheory = snrPerAntenna + bwGain + 10*log10(nantrx);
fprintf(1,'SNR theoretical = %7.2f\n', snrTheory);
%% Sweep the SNR
SNRtest = -10:5:50;
snrEq = zeros(size(SNRtest));
snrTheory = zeros(size(SNRtest));
nantrx = prod(nantUE);
bwGain = 10*log10(tx.waveformConfig.Nfft/tx.waveformConfig.NSubcarriers);
	
for idx = 1:length(SNRtest)
	snrPerAntenna = SNRtest(idx);
	Erx = mean(abs(y).^2, 'all');
	Enoise = 10.^(-0.1*snrPerAntenna)*Erx;
	ynoisy = y + sqrt(Enoise/2)*(randn(size(y)) +1i*randn(size(y)));

	rx.step(ynoisy);
	
	Err = mean(abs(rx.pdschSymRaw - rx.pdschChanEst.*tx.pdschSym).^2, 1);
	Erx = mean(abs(rx.pdschSymRaw).^2, 1);
	snrEq(idx) = mean(10*log10(Erx./Err));
	
	snrTheory(idx) = snrPerAntenna + bwGain + 10*log10(nantrx);
end

figure;
plot(SNRtest, snrEq, '-');
hold on;
plot(SNRtest, snrTheory, '-');
hold off;
grid on;
xlabel('SNR per Antenna');
ylabel('SNR');
legend('Equalized', 'Theory');
