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
phaseNoise = false;  % Include the phase noise of the RFFE
nonLin = true;		% Include the non-linearity of the RFFE
hpc = true;		% Run larger experiments using the NYU HPC 

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
%% ADC Calibration
% To set the noise power and ADC quantization levels, we run a large number
% of slots

if hpc
	ncal = 100;
else
	ncal = 10;
end
EyAvg = zeros(ncal, 1);

% Saturation levels of the lna to test
satLevTest = linspace(0, 20, 5);
numsat = length(satLevTest);

if ~hpc
	f = waitbar(0,'Calibrating');
end

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
    	
	for isat=1:numsat
		vsat = sqrt(mean(abs(y).^2, 'all')*10.^(0.1*satLevTest(isat)));
		lna = mmwsim.rffe.LNA('nonLin', nonLin, 'satLev', vsat);
		ymix = lna.step(y);
		
		% Measure average power
		EyAvg(i, isat) = mean(abs(ymix).^2, 'all');    
	end
	
	if ~hpc
		waitbar(i/ncal, f, 'Calibrating...');
	end
end
close(f);
EyAvg = mean(EyAvg);
%% Main simulation loop
% We next loop over the input SNRs and measure the SNR post equalization
% for each number of bitsrx
if hpc
	ADCtest = [4,6,0]';
	SNRtest = (-10:0.5:30)';
	satLevTest = 0:1:20;
	nit = 100;
else
	ADCtest = [4,6,0]';
	SNRtest = (-10:2:30)';
	satLevTest = linspace(0, 20, 5);
	nit = 10;
end

numsat = length(satLevTest);
numnbadc = length(ADCtest);
numsnr = length(SNRtest);
snrEq = zeros(numsnr, numsat, numnbadc, nit);

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
        
		for isat = 1:numsat
			% Add noise
			Enoise = 10.^(-0.1*snrAnt)*EyAvg(isat);
			ynoisy = y + sqrt(Enoise/2)*(randn(size(y)) +1i*randn(size(y)));

			for iadc = 1:numnbadc
				% Get the Saturation Level of the lna
				vsat = sqrt(mean(abs(y).^2, 'all')*10.^(0.1*satLevTest(isat)));
				
				% Get number of ADC bits
				nbadc = ADCtest(iadc);

				% Create the RX            
				rxInputVar = EyAvg(isat) + Enoise;
				rx = NRUERx(simParam, 'rxBF', wrx, 'nbitsADC', nbadc, ...
					'adcInputVar', rxInputVar, 'ncc', ncc, 'nonLin', nonLin, ...
					'phaseNoise', phaseNoise, 'satLev', vsat);
				rx.step(ynoisy);

				% Measure the post equalization SNR
				Err = mean(abs(rx.pdschSymRaw - rx.pdschChanEst.*tx.pdschSym).^2, 1);
				Erxeq = mean(abs(rx.pdschSymRaw).^2, 1);
				snrEq(isnr, isat, iadc, it) = mean(10*log10(Erxeq./Err));
			end
		end        
	end
    
	if ~hpc
		% Print progress
		for isat = 1:numsat
			for iadc = 1:numnbadc
				fprintf(1,'snr=%7.2f satLev=%7.2f nb=%d snrEq = %7.2f\n', ...
					snrAnt, satLevTest(isat), ADCtest(iadc), mean(snrEq(isnr,isat,iadc,:)));
			end
		end
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
	nplot = numsat + 1;
	legstr = cell(nplot,1);
	for i = 1:numsat
		satLev = satLevTest(i);
		legstr{i} = sprintf('satLev=%7.2f', satLev);
	end
	legstr{nplot} = 'Theory';

	for iadc = 1:numnbadc
		figure;
		nb = ADCtest(iadc);
		% Create plot
		snrEqAvg = mean(snrEq, 4);
		plot(SNRtest, snrEqAvg(:,:,iadc), '-', 'linewidth', 2);
		hold on;
		plot(SNRtest, snrTheory, '--', 'linewidth', 2);
		grid on;
		xlabel('SNR per Antenna');
		ylabel('SNR post-EQ');
		legend(legstr, 'Location', 'NorthWest');
		
		if nb == 0
			titlestr = sprintf('nbadc = inf | ncc = %d', ncc);
		else
			titlestr = sprintf('nbadc = %d | ncc = %d', nb, ncc);
		end
		
		if singPath
			titlestr = sprintf('%s | Single Path', titlestr);
		else
			titlestr = sprintf('%s | Multi Path', titlestr);
		end
		
		title(titlestr);
	end
end