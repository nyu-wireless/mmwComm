%% Test of AQN model for beamforming
% A simple program to estimate the effect of quantization and 
% compare this to a naive AQN model.
% 
% Test
% ----
% There are nx symbols and nrx antennas.  In the test, we:
%
% *  Generate a set of complex scalar symbols x(i), i=1,...,nx
% *  Pass them through a channel y(i,:) = x(i)*w(i,:) + d(i,:) where
%    w(i,:) is a random channel and d(i,:) is Gaussian noise.   
%    Note that y is nx x nrx
% *  Quantize this:   q = adc( y ).
% *  Estimate x from the quantized version:  xhat(i) = (q(i,:)*w')/(w'*w)
% *  Measure the error MSE = mean((x - xhat).^2)
% *  We plot the post-EQ SNR vs. SNR per antenna.
%       gam_in = SNR per ant = E|d(i,j)|^2
%       gam_out = post-EQ SNR = E|x(i)|^2/E|x(i) - xhat(i)|^2
%
% Native AQN model
% ----------------
% We compare the gam_out vs. gam_in curve to the prediction
% from a naive AQN model where we can treat the quantization noise as
% indepenent.
% * Suppose we can approximate the quantization error as:
%   q = adc(y) = y + dq,  dq = quantization error
% * Then q = x(i)*w(i,:) + d + dq
% * So, d + dq = effective noise
% * If we assume this effective noise is independent across antennas
%   we can predict the post-EQ SNR.
%
% What we see
% ------------
% If you set txSymType = 'iidPhase' and chanType = 'iidPhase', the 
% naive AQN model fits.
% If you set txSymType = 'iidGaussian' and chanType = 'iidPhase', the 
% naive AQN model does not fit.
% 
% Is there a more accurate model using the results in the paper.

%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
% We will use the following parameters
fc = 140e9;		% carrier frequency in Hz
nrx = 16;		% number of RX antennas
nx = 1e5;		% number samples per SNR point
nonLin = true;	% include the non-linearity of the RFFE
nit = 1;

snrTest = linspace(-10,40,10)';
nsnr = length(snrTest);

% ADC parameters 
nbits = 4;			% ADC bit-resolution
dither = true;		% randomize quantization error 
outType = 'float';
if fc == 28e9
	fsamp = 491.52e6;
elseif fc == 140e9
	fsamp = 1.9661e9;
end

% Transmit symbol type:  'iidGaussian' or 'iidPhase'
txSymType = 'iidGaussian';
xvar = 1;

% Channel type: 'iidGaussian', 'iidPhase', 'iidAoA' or 'ones'
chanType = 'iidGaussian';

% LNA parameters
lnaSatTest = linspace(5, 20, 4);	% Saturation level
lnaNoiseFig = 20;					% Noise figure
lnaGain = 1;						% Mixer gain
nlna = length(lnaSatTest);

% Mixer parameters
mixSatTest = 30;	% Saturation level
mixNoiseFig = 20;	% Noise figure
mixGain = 1;		% Mixer gain
nmix = length(mixSatTest);
%% Main simulation loop
% We next loop over the input SNRs and measure the SNR post equalization
% and we will compare this with the naive AQN model

snrEq = zeros(nsnr, nlna, nmix);
snrAntAQN = zeros(nsnr, nlna, nmix);
snrAQN = zeros(nsnr, nlna, nmix);

for it = 1:nit
	% Generate random channel
	if strcmp(chanType,'iidPhase')
		phase = 2*pi*rand(nx,nrx);
		w = exp(1i*phase);
	elseif strcmp(chanType, 'iidGaussian')
		w = (randn(nx,nrx) + 1i*randn(nx,nrx))/sqrt(2);
	elseif strcmp(chanType, 'randAoA')
		dsep = 0.5;
		theta = unifrnd(-pi/2,pi/2,nx,1);
		phase = 2*pi*cos(theta)*(0:nrx-1)*dsep;
		w = exp(1i*phase);
	elseif strcmp(chanType, 'ones')
		w = ones(nx,nrx);
	else
		error('Unknown channel type');
	end

	% Generate random symbols
	if strcmp(txSymType, 'iidGaussian')
		x = (randn(nx,1) + 1i*randn(nx,1))*sqrt(xvar/2);
	elseif strcmp(txSymType, 'iidPhase')
		x = exp(1i*2*pi*rand(nx,1));
	else
		error('Unknown TX symbol type');
	end

	% Generate RX symbols with no noise
	y0 = x.*w;

	% Loop over SNR values
	for isnr = 1:nsnr
		% Get the SNR and compute the noise variance
		snr = snrTest(isnr);
		dvar = xvar*10.^(-0.1*snr);

		% Add noise
		d = (randn(nx,nrx) + 1i*randn(nx,nrx))*sqrt(dvar/2);
		y = y0 + d;

		% Create and optimize the ADC
		yvar = xvar + dvar;
		adc = mmwsim.rffe.ADC('nbits', nbits, 'isComplex', true, ...
			'inputVar', yvar, 'dither', dither, 'outputType', outType);
		adc.optScale();            

		for ilna = 1:nlna
			% Get the LNA saturation level
			vsat = sqrt(mean(abs(y).^2, 'all')*10.^(0.1*lnaSatTest(ilna)));

			% LNA
			lna = mmwsim.rffe.LNA('nonLin', nonLin, 'satLev', vsat, ...
				'fsamp', fsamp,	'noiseFig', lnaNoiseFig);
			ylna = lna.step(y);

			% Get the mixer saturation level
			vsat = sqrt(mean(abs(y).^2, 'all')*10.^(0.1*mixSatTest));
			
			% Mixer
			mixer = mmwsim.rffe.Mixer('nonLin', nonLin, 'satLev', vsat, ...
				'fsamp', fsamp,	'noiseFig', mixNoiseFig);
			ymix = mixer.step(ylna);

			% Quantize
			q =  adc.step(ymix);
			adcVar = mean(abs(y-q).^2,'all');

			% Beamform
			xhat = sum(q.*conj(w),2) ./ sum(abs(w).^2,2);

			% Compute error
			xerr = mean(abs(x-xhat).^2);
			snrEq(isnr, ilna) = snrEq(isnr, ilna) + 10*log10(xvar / mean(abs(x-xhat).^2) );

			A = xvar^(-1)*(x'*xhat./nx).';
			T = (xhat' - A*x')*(xhat' - A*x')'./nx;

			% Naive AQN model
			snrAntAQN(isnr, ilna) = snrAntAQN(isnr, ilna) + 10*log10(xvar/(dvar + adcVar));        

			% AQN model
			snrAQN(isnr, ilna) = snrAQN(isnr, ilna) + 10*log10(abs(A)^2*xvar/T);
		end
	end
end
snrEq = snrEq./nit;
snrAntAQN = snrAntAQN./nit;
snrAQN = snrAQN./nit;

snrEqAQN = snrAntAQN + 10*log10(nrx);
%% Plot the results
for ilna=1:nlna
	figure;
	plot(snrTest, snrEq(:, ilna), 'o-', 'LineWidth', 2);
	hold on;
	plot(snrTest, snrEqAQN(:, ilna), '--', 'LineWidth', 2);
	plot(snrTest, snrAQN(:, ilna), '-x', 'LineWidth', 2);
	hold off;
	grid on;
	xlabel('SNR per antenna (dB)');
	ylabel('SNR post-EQ (dB)');
	legend('Simulation', 'Theory Naive AQN', 'Theory ISIT', 'Location', 'SouthEast');
	title(sprintf('fc = %d GHz, LNA: (sat = %d, NF = %d)\nnbadc = %d, Mixer: (sat = %d, NF = %d)', fc*1e-9, ...
		lnaSatTest(ilna), lnaNoiseFig, mixSatTest, mixNoiseFig, nbits));
end