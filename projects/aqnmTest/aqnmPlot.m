%% Load the simulation data
load('aqnmData1.mat');

%% Plot the Output SNR for all configurations

figure(1); clf;
hold on;
for isim=1:nsim
    plot(Pin, snrOut(:,isim), '-', 'linewidth', 1.5);
end
box on;
hold off;
axis tight;
xlabel('Receive power per antenna [dBm]', 'interpreter', 'latex', ...
    'fontsize', 13);
ylabel('Output SNR $\;(\gamma_\mathrm{out})\;$ [dB]', ...
    'interpreter', 'latex', ...
    'fontsize', 13);
%% Compare the Output SNR and the model

figure(1); clf;
cm = [0.5313, 0.7875, 0.7877; 0.3063 0.6166 0.6170; 0.1154 0.4519 0.4522];
Markers = {'o-', '-d', '-s', '-*', '->'};

theta = [100,300,600];
hold on;
for t=1:length(theta)
    plot(Pin, snrOut(:,theta(t)), 'color', cm(t,:), 'linewidth', 1.5);
    plot(Pin, rffeModel(:,theta(t)), Markers{t}, 'color', cm(t,:), ...
        'linewidth', 1.5, 'MarkerSize', 5);
end
box on;
hold off;
axis tight;
xlabel('Receive power per antenna [dBm]', 'interpreter', 'latex', ...
    'fontsize', 13);
ylabel('Output SNR $\;(\gamma_\mathrm{out})\;$ [dB]', ...
    'interpreter', 'latex', ...
    'fontsize', 13);
legend('$\mathrm{Sim}^{(1)}$', '$\mathrm{Model}^{(1)}$', ...
    '$\mathrm{Sim}^{(2)}$','$\mathrm{Model}^{(2)}$',...
    '$\mathrm{Sim}^{(3)}$', '$\mathrm{Model}^{(3)}$',...
    'interpreter', 'latex', 'location', 'best', 'Box','off', ...
    'fontsize', 11);

%% Find the minimum Power for each configuration

nbins = 200;
[~, edge0, ~] = histcounts(rffeNF, nbins);
[~, edge1, ~] = histcounts(snrSat, nbins);

p = zeros(nbins-1, nbins-1);
figure;
for i=1:nbins-1
    for j=1:nbins-1
        I = (rffeNF <= edge0(i+1)) & (snrSat >= edge1(j));
        if (sum(I) > 0)
            p(j, nbins-i) = min(rffePower(I),[],'all');
        else
            p(j, nbins-i) = p(j-1, nbins-i);
        end
    end
end
imagesc(rffeNF, snrSat, p)

ylabel('Target saturation SNR $\;(\gamma_\mathrm{sat}^\mathrm{tgt})\;$ [dB]', ...
    'interpreter', 'latex', ...
    'fontsize', 13);
xlabel('Target noise figure $\;(F^\mathrm{tgt})\;$ [dB]', ...
    'interpreter', 'latex', 'fontsize', 13);
c = colorbar;
set(c, 'fontsize', 13);
c.Label.String = 'Power Consumption [mW]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'YDir','normal');

%% Plot the NF, Gain and IIP3 of the IF Mixers
cm = [0.53, 0.78, 0.78; 0.42, 0.70, 0.70;
    0.31, 0.62, 0.62; 0.21, 0.58, 0.58; 0.12, 0.45, 0.45];

Markers = {'o-', '-d', '-s', '-*', '->'};

% Noise Figure
figure(1);clf;
hold on;
for imix=1:nmix
    plot(mixPLO, mixNF(:,imix), Markers{imix}, 'color', cm(imix,:), ...
        'linewidth', 1.5, 'markersize', 5);
end
hold off;
ylabel('Noise Figure [dB]', 'interpreter', 'latex', ...
    'fontsize', 17);
xlabel('Power of the Local Oscillator [dBm]', 'interpreter', 'latex', ...
    'fontsize', 17);
l=legend(...
    '$\mathrm{Mixer}^{(1)}$', '$\mathrm{Mixer}^{(2)}$', '$\mathrm{Mixer}^{(3)}$',...
    '$\mathrm{Mixer}^{(4)}$', '$\mathrm{Mixer}^{(5)}$', ...
    'location', 'northoutside', 'orientation', 'horizontal',...
    'numcolumns', 4, 'Box','off', 'fontsize', 14);
set(l, 'interpreter', 'latex')
box on

% Gain
figure(2);clf;
hold on;
for imix=1:nmix
    plot(mixPLO, mixGain(:,imix), Markers{imix}, 'color', cm(imix,:), ...
        'linewidth', 1.5, 'markersize', 5);
end
hold off;
ylabel('Gain [dB]', 'interpreter', 'latex', ...
    'fontsize', 17);
xlabel('Power of the Local Oscillator [dBm]', 'interpreter', 'latex', ...
    'fontsize', 17);
l=legend(...
    '$\mathrm{Mixer}^{(1)}$', '$\mathrm{Mixer}^{(2)}$', '$\mathrm{Mixer}^{(3)}$',...
    '$\mathrm{Mixer}^{(4)}$', '$\mathrm{Mixer}^{(5)}$', ...
    'location', 'northoutside', 'orientation', 'horizontal',...
    'numcolumns', 4, 'Box','off', 'fontsize', 14);
set(l, 'interpreter', 'latex')
box on;

% IIP3
figure(3);clf;
hold on;
for imix=1:nmix
    plot(mixPLO, mixIIP3(:,imix), Markers{imix}, 'color', cm(imix,:), ...
        'linewidth', 1.5, 'markersize', 5);
end
hold off;
ylabel('IIP3 [dBm]', 'interpreter', 'latex', ...
    'fontsize', 17);
xlabel('Power of the Local Oscillator [dBm]', 'interpreter', 'latex', ...
    'fontsize', 17);
l=legend(...
    '$\mathrm{Mixer}^{(1)}$', '$\mathrm{Mixer}^{(2)}$', '$\mathrm{Mixer}^{(3)}$',...
    '$\mathrm{Mixer}^{(4)}$', '$\mathrm{Mixer}^{(5)}$', ...
    'location', 'northoutside', 'orientation', 'horizontal',...
    'numcolumns', 4, 'Box','off', 'fontsize', 14);
set(l, 'interpreter', 'latex')
box on;