%% 5G NR downlink PDSCH - Non Linear RFFE

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
snrInTest = linspace(-10, 50, 21)';

% ADC resolution (4-bit, 5-bit, 6-bit). For inf-bit use 0
adcTest = [4, 5, 6];

% Define a uniform random function in [a,b]
rnd = @ (a,b) (a + (b-a)*rand(1,1));

% Load the RFFE models
if fc == 140e9
	load('rffe_140GHz.mat');
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
	nit = 20;	% number of iterations
else
	rng('shuffle');
	arrayId = 0;
	nit = 1;	% number of iterations
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
ndrivers = length(0:log2(nrx));		% num of LO drivers
nsim = nlna * nmix * nplo * nadc;	% num of simulations

% Intialize the output SNR vector
snrOut = zeros(nsnr, nsim, nit);

% Loop over the slots
for it = 1:nit
tic;
	fprintf('\nSlot: %d\n',it);
	
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
	
	% Rescale so that it is Es/kT = 1
	scale = sqrt(EkT/mean(abs(y).^2, 'all'));
	y = y * scale;
			
	% Create all possible receiver configurations
	isim = 1;
	rx = cell(nsim, 1);
	for iadc = 1:nadc
		for ilna = 1:nlna
			for imix = 1:nmix
				for iplo = 1:nplo
					rx{isim} = NRUERx(...
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
					isim = isim + 1;
				end
			end
		end
	end
	
	% If it is possible run parrallel simulations to calculate the output
	% SNR of a linear receiver
	parfor isim = 1:nsim
		snrOut(:, isim, it) = rx{isim}.step(y);
	end
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
	rffePower(isim, :) = rx{isim}.power();	% RFFE power consumption [mW]
	rffeNF(isim) = rx{isim}.nf();			% Effective noise figure [dBm]
end

% Find the minimum power for each theta based on the different LO drivers.
[rffePower, idriver] = min(rffePower, [], 2);

%% Fit a model
% We use two parameters to characterize the performance of each
% configuration: (a) effective noise figure that is dominant in low input
% power; (b) saturation SNR that is dominant in high SNR. Using these
% values we can fit a model for the output SNR as follows,

% For non-linear system we empirically show that the output SNR can be
% calculated by the following formula
nom = nrx * 10.^(0.1*snrInTest);
denom = reshape(10.^(0.1*rffeNF), [], nsim) + ...
	nrx * reshape(10.^(0.1*snrInTest), nsnr, []) .* ...
	reshape(10.^(-0.1*snrSat), [], nsim);
rffeModel = 10*log10(nom./denom);

%% Save the data

% clear the objects and some of the variables
clear tx chan rx simParam arrgNB arrgNB0 arrUE arrUE0 elem elemInterp rnd ...
	x y aoaEl aoaAz aodAz aodEl it isim isnr ilna imix iplo iadc nom ...
	denom noiseTemp EkT lambda len scale vp elemGainRx elemGainTx dly ...
	dsep gain urx wrx utx wtx groundPlaneLen dlySpread;

% save the simulation results
save(sprintf('nrData_%d.mat', arrayId));