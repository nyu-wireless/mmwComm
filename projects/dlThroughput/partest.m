%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
fc = 140e9;  % carrier frequency in Hz
noiseTemp = 290;  % noise temperature in K
EkT = physconst('Boltzman')*noiseTemp;  % noise energy
isLOS = true;  % true to select a line-of-slight channel with no multipath
isFD = false;  % true to select fully-digital. False for analog beamforming.
isLinear = false;  % true to select a linear RF front-end

if fc == 140e9
    load('../../models/rffe140GHz.mat')
else
    load('../../models/rffe28GHz.mat')
end

ueVel = [5 0 0];	% UE velocity vector in m/s
nantUE = [4,4];		% array size at the UE (mobile device)
nantgNB = [8,8];	% array size at the gNB (base station)
nrx = nantUE(1)*nantUE(2);	% num of receive antenna

if fc == 140e9
    SubcarrierSpacing = 240;  % SCS in KHz
else
    SubcarrierSpacing = 120;  % SCS in KHz
end
NRB = 66;  % num of resource blocks
nscPerRB = 12;  % num of sub-carrier per resource block
Modulation = 'QPSK'; % modulation, i.e., 'QPSK', '16QAM', '64QAM'
ncc = 1;  % num of component carriers

% Define a uniform random function in [a,b]
rnd = @ (a,b) (a + (b-a)*rand(1,1));

%% Create the 3GPP 5G NR objects
carrierConfig = nrCarrierConfig(...
    'NSizeGrid', NRB, ...
    'SubcarrierSpacing', SubcarrierSpacing);

% Create the waveform configuration from the carrier configuration
waveformConfig = nrOFDMInfo(carrierConfig);

% Scale the bandwidth based on the number of component carriers
waveformConfig.SampleRate = ncc * waveformConfig.SampleRate;
waveformConfig.ncc = ncc;
waveformConfig.fcc = [0, 0];

ptrsConfig = nrPDSCHPTRSConfig(...
    'TimeDensity', 1, ...
    'FrequencyDensity', 2, ...
    'REOffset', '00');

pdschConfig = nrPDSCHConfig(...
    'Modulation', Modulation, ...
    'PRBSet', 0:NRB-1, ...
    'SymbolAllocation', [1, carrierConfig.SymbolsPerSlot-1], ...
    'EnablePTRS', 1, ...
    'PTRS', ptrsConfig);

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

%% Create the RFFE configurations
adcTest = [4,5,6,7,8];
idsn = 1;
ndsn = length(lna) * size(mixer,1) * size(mixer, 2) * ...
    length(adcTest) * (1 + length(ps));

rffe = cell(ndsn,1);
for ilna = 1:length(lna)
    for imix = 1:size(mixer, 1)
        for iplo = 1:size(mixer, 2)
            for iadc = 1:length(adcTest)
                rffe{idsn} = {lna{ilna}, mixer{imix, iplo}, adcTest(iadc), {}};
                idsn = idsn + 1;
                for ips = 1:length(ps)
                    rffe{idsn} = {lna{ilna}, mixer{imix, iplo}, adcTest(iadc), ps{ips}};
                    idsn = idsn + 1;
                end
            end
        end
    end
end
%%
clc
nsig = 13;
nslot = 5;
nf = zeros(ndsn, 1);
power = zeros(ndsn, 1);
idx = zeros(ndsn, 1);
snrOut = zeros(nsig, ndsn, nslot);

for islot = 1:nslot
    if isLOS
        gain = 0;
        aoaAz = rnd(-180, 180);
        aoaEl = rnd(-90, 90);
        aodAz = rnd(-180, 180);
        aodEl = rnd(-90, 90);
        dly = 0;
        im = 1;
    end

    [utx, ~] = arrgNB.step(aodAz, aodEl);
    [urx, ~] = arrUE.step(aoaAz, aoaEl);


    wtx = conj(utx(:,im)); % TX direction at the gNB
    wtx = wtx / norm(wtx);

    wrx = conj(urx(:,im)); % RX direction at the UE
    wrx = wrx / norm(wrx);


    gNB = NRgNBTx(...
        'carrierConfig', carrierConfig, ...
        'pdschConfig', pdschConfig, ...
        'waveformConfig', waveformConfig);

    % Create the channel
    chan = mmwsim.chan.MIMOMPChan(...
        'aoaAz', aoaAz, ...
        'aoaEl', aoaEl, ...
        'aodAz', aodAz, ...
        'aodEl', aodEl, ...
        'gain', gain, ...
        'rxArr', arrUE, ...
        'txArr', arrgNB, ...
        'dly', dly, ...
        'fsamp', waveformConfig.SampleRate);

    % Create TX and RX
    z = gNB.step();
    z = z*wtx';  % Apply TX beamforming
    z = z - mean(z, 'all');  % Make sure the TX signal is zero mean
    y = chan.step(z);  % pass the TX signal through a channel

    % Rescale so that it is Es/(fs*kT) = 1
    scale = sqrt((waveformConfig.SampleRate*EkT)/mean(abs(y).^2, 'all'));
    y = y * scale;
    x = gNB.pdschSym;  % TX known signal to measure SNR
    tic
    parfor idsn = 1:ndsn
        lna0 = rffe{idsn}{1};
        mix0 = rffe{idsn}{2};
        adc0 = rffe{idsn}{3};
        ps0 = rffe{idsn}{4};
        nstage0 = rffe{idsn}{5};
        isFD = isempty(ps0);

        ue = NRUERx(...
            'carrierConfig', carrierConfig, ...
            'pdschConfig', pdschConfig, ...
            'waveformConfig', waveformConfig, ...
            'ps', ps0, ...
            'lna', lna0, ...
            'nstages', nstage0, ...
            'mixer', mix0, ...
            'nbits', adc0, ...
            'isLinear', isLinear, ...
            'isFD', isFD, ...
            'nrx', nrx);

        snrSigTest = linspace(-10, 50, nsig);
        for isig = 1:nsig
            % Scale the signal power
            xnoisy = y .* 10^(0.05*snrSigTest(isig));

            % The output SNR at that signal level
            rxPDSCHSym = ue.step(xnoisy, wrx);
            xhat = rxPDSCHSym;
            xvar = mean(abs(x).^2, 'all');
            a = mean(conj(xhat).*x)/mean(abs(x).^2);
            dvar = mean(abs(xhat - a*x).^2);
            snrOut(isig, idsn, islot) = 10*log10(abs(a).^2*xvar/dvar);
        end

        ue.rffe.SampleRate = waveformConfig.SampleRate*8;
        nf(idsn) = ue.rffe.nf();
        [power(idsn), idx(idsn)] = min(ue.rffe.power());
    end
    toc
end

%%
snrSigTest = linspace(-10, 50, nsig);
snrOut0 = mean(snrOut, 3);

figure(2);
plot(snrSigTest, snrOut0, '-o');
axis tight;
xlabel('Input SNR $\;(\gamma_\mathrm{in})\;$ [dB]', 'interpreter', 'latex', ...
    'fontsize', 13);
ylabel('Output SNR $\;(\gamma_\mathrm{out})\;$ [dB]', ...
    'interpreter', 'latex', ...
    'fontsize', 13);
legend('Fully-Digital', 'Analog (Active)', ...
    'Analog (Passive -15 dB)', 'Analog (Passive -10 dB)', 'location', 'best');

for idsn = 1:ndsn
    fprintf(1, "Pdc: %.0f mW\tNF: %.2f dB\n", power(idsn), nf(idsn));
end
%%
snrSat = max(snrOut0);
nom = nrx * 10.^(0.1*snrSigTest).';
denom = reshape(10.^(0.1*nf), [], ndsn) + ...
    nrx * reshape(10.^(0.1*snrSigTest), nsig, []) .* ...
    reshape(10.^(-0.1*snrSat), [], ndsn);
snrSat = snrSat.';
rffeModel = 10*log10(nom./denom);
%%
k = 1;
plot(snrSigTest, snrOut(:,k));
hold on;
plot(snrSigTest, rffeModel(:,k));
hold off;

%% Find the minimum Power for each configuration

nbins = 200;
[~, edge0, ~] = histcounts(nf, nbins);
[~, edge1, ~] = histcounts(snrSat, nbins);

p = zeros(nbins-1, nbins-1);
for i=1:nbins-1
    for j=1:nbins-1
        I = (nf <= edge0(i+1)) & (snrSat >= edge1(j));
        if (sum(I) > 0)
            p(j, nbins-i) = min(power(I),[],'all');
        else
            p(j, nbins-i) = p(j-1, nbins-i);
        end
    end
end

%%
[X, Y] = meshgrid(edge0(end:-1:3), edge1(1:end-2));
s = surf(X, Y, p, 'EdgeColor', [0 0 0], 'LineStyle', '-', 'FaceLighting', 'gouraud',...
    'FaceAlpha', 1, 'EdgeAlpha', 0.2);

ylabel('Target saturation SNR $\;(\gamma_\mathrm{sat}^\mathrm{tgt})\;$ [dB]', ...
    'interpreter', 'latex', ...
    'fontsize', 13);
xlabel('Target noise figure $\;(F^\mathrm{tgt})\;$ [dB]', ...
    'interpreter', 'latex', 'fontsize', 13);
zlabel('Power Consumption [mW]', ...
    'interpreter', 'latex', 'fontsize', 13);
c = colorbar;
set(c, 'fontsize', 13);
c.Label.String = 'Power Consumption [mW]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'YDir','normal');
xh = get(gca, 'XLabel');
set(xh, 'Units', 'Normalized');
pos = get(xh, 'Position');
set(xh, 'Position', pos.*[1.05, -0.5, 1], 'Rotation', 14);