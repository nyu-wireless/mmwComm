%hCarrierAggregationPlotSpectrum Plot spectrum of waveform
%   hCarrierAggregationPlotSpectrum(WAVE,SR,TITLE,LEGEND) plots the
%   spectrum of a waveform WAVE at sampling rate SR using the spectrum
%   analyzer. TITLE and LEGEND, specify the title and legend of the plot.
%   For multiple plots, LEGEND can be a cell array of character vectors.

%   Copyright 2015-2016 The MathWorks, Inc.

function specPlot = hCarrierAggregationPlotSpectrum(waveform,SR,title,legend)

specPlot = dsp.SpectrumAnalyzer;
specPlot.SampleRate = SR;
% The window length controls the frequency resolution 
specPlot.FrequencyResolutionMethod = 'WindowLength';
% FFT length expressed as a property in FFTLength
specPlot.FFTLengthSource = 'Property';
fftLength = (SR/30.72e6*2048);
specPlot.FFTLength = fftLength;
% calculate number of spectral averages in spectrum
numSpectralAverages = fix(size(waveform, 1)/fftLength);
specPlot.SpectralAverages = numSpectralAverages;
% overlap between spectral averages
specPlot.OverlapPercent = 0;
% Set windowing parameters
% specPlot.Window = 'Rectangular';
specPlot.WindowLength = fftLength;
specPlot.Window = 'Kaiser';
% set legend and title
specPlot.ChannelNames = legend;
specPlot.ShowLegend = true;
specPlot.Title = title;
specPlot.YLimits = [-120 0];
% plot spectrum
specPlot(waveform(1:(fftLength*numSpectralAverages),:));

end

