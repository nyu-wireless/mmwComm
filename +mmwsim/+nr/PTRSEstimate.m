function [ofdmGrid] = ptrsEstimate(ptrsIndices, ptrsSym, ofdmGrid)
	nptrs = 33;	% Number of PTRS per OFDM Symbol
	I = unique(floor(ptrsIndices/792))+1;
	
	% phaseDiff = angle(mean(reshape(ofdmGrid(ptrsIndices).*conj(ptrsSym), nptrs, [])));
	phaseDiff = mean(unwrap(angle(reshape(ofdmGrid(ptrsIndices).*conj(ptrsSym), nptrs, []))));
	phase = interp1(I, phaseDiff, 1:14, 'pchip');
	
	ofdmGrid = ofdmGrid .* exp(-1j*phase);
end