clear; clc; setup; config_irs;

%% ! R-E region for fixed and adaptive IRS
reAdaptiveIrsSample = cell(nChannels, length(Variable.bandwidth));
reFsIrsSample = cell(nChannels, length(Variable.bandwidth));
reWitIrsSample = cell(nChannels, length(Variable.bandwidth));
reWptIrsSample = cell(nChannels, length(Variable.bandwidth));
reNoIrsSample = cell(nChannels, length(Variable.bandwidth));

reAdaptiveIrsSolution = cell(nChannels, length(Variable.bandwidth));
reFsIrsSolution = cell(nChannels, length(Variable.bandwidth));
reWitIrsSolution = cell(nChannels, length(Variable.bandwidth));
reWptIrsSolution = cell(nChannels, length(Variable.bandwidth));
reNoIrsSolution = cell(nChannels, length(Variable.bandwidth));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

	for iBandwidth = 1 : length(Variable.bandwidth)
		% * Calculate carrier frequency
		[subbandFrequency] = subband_frequency(centerFrequency, Variable.bandwidth(iBandwidth), nSubbands);

		% * Construct channels
		[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
		[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
		[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

		% * Adaptive IRS and waveform design
		[reAdaptiveIrsSample{iChannel, iBandwidth}, reAdaptiveIrsSolution{iChannel, iBandwidth}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		witCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{1}.compositeChannel;
		wptCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{end}.compositeChannel;

		% * Upper bound by frequency-selective IRS
		[fsIrsCompositeChannel] = fs_irs_composite_channel(directChannel, incidentChannel, reflectiveChannel);
		[reFsIrsSample{iChannel, iBandwidth}, reFsIrsSolution{iChannel, iBandwidth}] = re_sample_reference(beta2, beta4, fsIrsCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization with fixed WIT-optimized IRS
		[reWitIrsSample{iChannel, iBandwidth}, reWitIrsSolution{iChannel, iBandwidth}] = re_sample_reference(beta2, beta4, witCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization with fixed WPT-optimized IRS
		[reWptIrsSample{iChannel, iBandwidth}, reWptIrsSolution{iChannel, iBandwidth}] = re_sample_reference(beta2, beta4, wptCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization without IRS
		[reNoIrsSample{iChannel, iBandwidth}, reNoIrsSolution{iChannel, iBandwidth}] = re_sample_reference(beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
	end
end

% * Average over channel realizations
reAdaptiveIrsInstance = cell(1, length(Variable.bandwidth));
reFsIrsInstance = cell(1, length(Variable.bandwidth));
reWitIrsInstance = cell(1, length(Variable.bandwidth));
reWptIrsInstance = cell(1, length(Variable.bandwidth));
reNoIrsInstance = cell(1, length(Variable.bandwidth));

for iBandwidth = 1 : length(Variable.bandwidth)
	reAdaptiveIrsInstance{iBandwidth} = mean(cat(3, reAdaptiveIrsSample{:, iBandwidth}), 3);
	reFsIrsInstance{iBandwidth} = mean(cat(3, reFsIrsSample{:, iBandwidth}), 3);
	reWitIrsInstance{iBandwidth} = mean(cat(3, reWitIrsSample{:, iBandwidth}), 3);
	reWptIrsInstance{iBandwidth} = mean(cat(3, reWptIrsSample{:, iBandwidth}), 3);
	reNoIrsInstance{iBandwidth} = mean(cat(3, reNoIrsSample{:, iBandwidth}), 3);
end

% * Save batch data
save(sprintf('data/re_irs/re_irs_%d.mat', iBatch));
