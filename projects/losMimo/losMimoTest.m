%% LOS MIMO capacity simulation
% The program estimates the distribution of rates using high-rank LOS MIMO
% on UAV-gNB link.  The simulation is as follows:
%
% * A gNB ground cell is fixed at the origin with an upward pointing
%   antenna
% * A UAV (UE) is randomly placed above the ground cell
% * We compute the channel from the UAV to the gNB (uplink)
% * Compute the rate of the channel
% * The CDF of the rate is plotted over multiple trials
%
% In the future, we want to add channel estimation error

%% Packages
% Add the folder containing +mmwsim to the MATLAB path.
addpath('../..');

%% Parameters
% We will use the following parameters
fc = 140e9;			% carrier frequency in Hz
nantUE = [4,4];		% array size at the UE (mobile device)
nantgNB = [8,8];	% array size at the gNB (base station)
arrAppUE = [0.4 0.4];    % array aperture at the UE (m)
arrAppgNB = [0.4 0.4];   % array aperture at the gNB (m)

% Position of the UAV is randomly placed in this box
uePosMin = [-50,-50,50];
uePosMax = [50,50,150];

EsN0 = -10;       % SNR per sample per antenna in dB
nstreamsMax = 4;  % maximum number of streams

nit = 100;  % number of trials per distance

% Link budget parameters
kT = -174;  % Thermal noise kT
Ptx = 26;   % Transmit power
fsamp = 1.6e9;  % Sample rate in Hz
noiseFig = 6;     % noise figure


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
dsepUE = arrAppUE ./ (nantUE-1);
dsepgNB = arrAppgNB ./ (nantgNB-1);
arrUE0  = phased.URA(nantUE,dsepUE,'ArrayNormal','x');
arrgNB0 = phased.URA(nantgNB,dsepgNB,'ArrayNormal','x');


% Create a wrapper class around the arrays to handle orientation modeling.
% The ArrayPlatform class includes an array along with axes to store the
% orientation of the array relative to some global coordinate system.
arrgNB = mmwsim.antenna.ArrayPlatform('arr', arrgNB0, 'elem', elemInterp, ...
    'fc', fc);
arrUE  = mmwsim.antenna.ArrayPlatform('arr', arrUE0, 'elem', elemInterp, ...
    'fc', fc);

% Align the arrays so that the gNB is pointing directly upwards
% and the UE is pointing downwards
arrgNB.alignAxes(0,90);
arrUE.alignAxes(0,-90);

%% Create a MIMO rate capacity object
% This object is used to estimate the capacity at each point
rateCalc = mmwsim.rate.MIMORateCalc('includeRateLoss', true, ...
    'nstreamsMax', nstreamsMax);


%% Compute the TX SNR
% We next compute Etx/N0, the total TX energy per symbol divided
% by the noise energy per sample
EtxN0 = Ptx - 10*log10(fsamp) - kT - noiseFig;


%% Main simulation loop

% Generate random UAV positions
uePos = rand(nit,3).*(uePosMax-uePosMin) + uePosMin;

rate = zeros(nit,1);
for it = 1:nit
    
    % Generate a random position of the UE
    arrUE.set('pos', uePos(it,:) );
    
    % Find the narrowband channel
    chan = mmwsim.chan.LOSMIMOChan('rxArr', arrUE, 'txArr', arrgNB, ...
        'fc', fc);
    chan.computeChanMatrix();
    H = chan.chanMatrix;
    
    % Compute the rate
    rate(it) = rateCalc.computeRate(H,EtxN0);
end

%% Plot the CDF of the rate
p = (1:nit)/nit;
plot(sort(rate),p, 'linewidth', 3);
grid on;
xlabel('Rate (bits/samp)');
ylabel('CDF');