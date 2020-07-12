%% Test of AQN model for beamforming
% A simple program to estimate the effect of quantization.
% There are nx symbols and nrx antennas
%
% *  Generate a set of symbols x(i), i=1,...,nx
% *  Pass them through a channel y(i,:) = x(i)*w(i,:) + d(i,:) where
%    w(i,:) is a random channel.   Note that y is nx x nrx
% *  Quantize this:   q = adc( y ).  
% *  Estimate x from the quantized version:  xhat(i) = (q(i,:)*w')/(w'*w)
% *  Measure the error MSE = mean((x - xhat).^2)
%
% Simple theory:
% * Suppose we can approximate the quantization error as:
%   q = adc(y) = y + dq,  dq = quantization error
% * Then q = x(i)*w(i,:) + d + dq
% * So, d + dq = effective noise
% * If we assume this effective noise is independent across antennas
%   we can predict the post-EQ SNR. 
% But, this prediction does not work.  

nbits = 4;
nrx = 16;   % number of RX antennas
nx = 10000;  % number samples per SNR point
snrTest = linspace(0,50,10)';
nsnr = length(snrTest);
xvar = 1;
phaseDither = false;

snrEq = zeros(nsnr,1);
snrAntAQN = zeros(nsnr,1);

for isnr = 1:nsnr
    % Generate random symbols
    x = (randn(nx,1) + 1i*randn(nx,1))*sqrt(xvar/2);

    % Get the SNR and compute the noise variance
    snr = snrTest(isnr);
    dvar = xvar*10.^(-0.1*snr);

    % Create and optimize the ADC
    yvar = xvar + dvar;
    adc = mmwsim.rffe.ADC('nbits', nbits, 'isComplex', true, ...
        'inputVar', yvar, 'phaseDither', phaseDither);
    adc.optScale();

    % Random channel
    dsep = 0.5;
    theta = unifrnd(-pi/2,pi/2,nx,1);
    phase = 2*pi*cos(theta)*(0:nrx-1)*dsep;
    %phase = 2*pi*rand(nx,nrx);
    w = exp(1i*phase);
    y0 = x.*w;

    % Add noise    
    d = (randn(nx,nrx) + 1i*randn(nx,nrx))*sqrt(dvar/2);
    y = y0 + d;

    % Quantize
    q = adc.step(y);

    % Beamform
    xhat = sum(q.*conj(w),2) ./ sum(abs(w).^2,2);
    xerr = mean(abs(x-xhat).^2);
    snrEq(isnr) = 10*log10(xvar / mean(abs(x-xhat).^2) );

    % AQN model
    snrAntAQN(isnr) = 10*log10(xvar/(dvar + adc.quantVar));
        
end

snrEqAQN = snrAntAQN + 10*log10(nrx);
plot(snrTest, snrEq, 'o-');
hold on;
plot(snrTest, snrEqAQN, '--');
hold off;
grid on;

% AQN model



