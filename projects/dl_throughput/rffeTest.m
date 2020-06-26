%% Test of the RFFE
% Test of the RFFE performance using on the 5G NR downlink PDSCH.
% The post equalized SNR is measured as a function of the input SNR
% and the RFFE parameters.  Right now, we only test the input ADC
% resolution.

%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
% We will use the following parameters
fc = 28e9;          % carrier frequency in Hz
nantUE = [4,4];     % array size at the UE (mobile device)
nantgNB = [8,8];    % array size at the gNB (base station)
ueVel = [5 0 0];    % UE velocity vector in m/s
ncc = 1;            % number of component carriers (1 or 4)
singPath = true;	% Use only the strongest path
phaseNoise = false; % Include the phase noise of the RFFE
nonLin = true;		% Include the non-linearity of the RFFE
hpc = false;		% Run larger experiments using the NYU HPC 

% Parameters describing the PDSCH  
simParam = PDSCHSimParam('fc', fc);

%% NYU HPC
if hpc
	aID = getenv('SLURM_ARRAY_TASK_ID');
	if isempty(aID)
		aID = '1';
	end
	rng(str2double(aID),'twister');
else
	rng('shuffle');
end

%% Load the 3GPP NR channel model
% We load the channel using the MATLAB toolbox
if singPath
    gain = 0;
    aoaAz = 0;
    aoaEl = 0;
    aodAz = 0;
    aodEl = 0;
    dly = 0;
else
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
end

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

% Create wrapper class
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

% Find the paths with maximum gain
[gainmax, im] = max(gain);

% Assume the arrays are aligned on the maximum path
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

%% Calibration
% To set the noise power and ADC quantization levels, we run a large number
% of slots
ncal = 10;
EyAvg = zeros(ncal,1);

f = waitbar(0,'Calibrating');
for i = 1:ncal
    
    % Create a TX object using the NRgNBTx object for each CC.
    tx = NRgNBTx(simParam, 'nonLin', nonLin, 'phaseNoise', phaseNoise);
    tx.set('txBF', wtx, 'ncc', ncc);
    
    % Create the channel
    chan = mmwsim.chan.MIMOMPChan('aoaAz', aoaAz, 'aoaEl', aoaEl, 'aodAz', aodAz, ...
        'aodEl', aodEl, 'gain', gain, 'rxArr', arrUE, 'txArr', arrgNB, ...
        'dly', dly, 'fsamp', tx.waveformConfig.SamplingRate);
    
    % Simulate TX and channel
    x = tx.step();
    y = chan.step(x);
    
    % Measure average power
    EyAvg(i) = mean(abs(y).^2, 'all');    
    
    waitbar(i/ncal,f, 'Calibrating...');
    
end
close(f);
EyAvg = mean(EyAvg);

%% Main simulation loop
% We next loop over the input SNRs and measure the SNR post equalization
% for each number of bitsrx
if hpc
	ADCtest = [0,4,6]';
	SNRtest = (-10:0.5:30)';
	nit = 100;
else
	ADCtest = [0,4,6]';
	SNRtest = (-10:2:30)';
	nit = 10;
end

numnbadc = length(ADCtest);
numsnr = length(SNRtest);
snrEq = zeros(numsnr, numnbadc, nit);

% Loop over SNR values
for isnr = 1:numsnr
    snrAnt = SNRtest(isnr);
    
    % Loop over the slots
    for it=1:nit
        
        % Create a TX object using the NRgNBTx object for each CC.
        tx = NRgNBTx(simParam, 'nonLin', nonLin, 'phaseNoise', phaseNoise);
        tx.set('txBF', wtx, 'ncc', ncc);
        
        % Create the channel
        chan = mmwsim.chan.MIMOMPChan('aoaAz', aoaAz, 'aoaEl', aoaEl, ...
			'aodAz', aodAz, 'aodEl', aodEl, 'gain', gain, 'rxArr', arrUE, ...
			'txArr', arrgNB, 'dly', dly, 'fsamp', tx.waveformConfig.SamplingRate);
        
        % Create TX and RX
        x = tx.step();
        y = chan.step(x);
        
        % Add noisex
        Enoise = 10.^(-0.1*snrAnt)*EyAvg;
        ynoisy = y + sqrt(Enoise/2)*(randn(size(y)) +1i*randn(size(y)));
        
        for iadc = 1:numnbadc
            
            % Get number of ADC bits
            nbadc = ADCtest(iadc);
            
            % Create the RX            
            rxInputVar = EyAvg + Enoise;
            rx = NRUERx(simParam, 'rxBF', wrx, 'nbitsADC', nbadc, ...
                'adcInputVar', rxInputVar, 'ncc', ncc, 'nonLin', nonLin, ...
				'phaseNoise', phaseNoise);      
            rx.step(ynoisy);
            
            % Measure the post equalization SNR
            Err = mean(abs(rx.pdschSymRaw - rx.pdschChanEst.*tx.pdschSym).^2, 1);
            Erxeq = mean(abs(rx.pdschSymRaw).^2, 1);
            snrEq(isnr,iadc,it) = mean(10*log10(Erxeq./Err));
        end
        
    end
    
    % Print progress
    for iadc = 1:numnbadc
        fprintf(1,'snr=%7.2f nb=%d snrEq = %7.2f\n', snrAnt, ADCtest(iadc),...
            mean(snrEq(isnr,iadc,:)) );
    end
end

%% Compute theoretical maximum SNR
nantrx = prod(nantUE);
bwGain = 10*log10(tx.waveformConfig.Nfft/tx.waveformConfig.NSubcarriers);
snrTheory = SNRtest + bwGain + 10*log10(nantrx);

%% Plot the results
if hpc
	% Save the workspace instead of ploting
	save(strcat('workspace','_',num2str(aID),'.mat'))
else
	% Create legend strings
	nplot = numnbadc + 1;
	legstr = cell(nplot,1);
	for i = 1:numnbadc
		nb = ADCtest(i);
		if nb == 0
			legstr{i} = 'nbits=inf';        
		else
			legstr{i} = sprintf('nbits=%d', nb);        
		end
	end
	legstr{nplot} = 'Theory';

	% Create plot
	snrEqAvg = mean(snrEq,3);
	plot(SNRtest, snrEqAvg, '-', 'linewidth', 2);
	hold on;
	plot(SNRtest, snrTheory, '--', 'linewidth', 2);
	grid on;
	xlabel('SNR per Antenna');
	ylabel('SNR post-EQ');
	legend(legstr, 'Location', 'NorthWest');
	title('Single-Carrier / Single-Path');
end