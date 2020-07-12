function out = hCarrierAggregationModulate(in, sr, f)

    % Create vector of time samples
    t = ((0:size(in,1)-1)/sr).';   
    
    % For each antenna, apply the frequency offset correction
	out = in .* exp(1i*2*pi*f*t);
end
