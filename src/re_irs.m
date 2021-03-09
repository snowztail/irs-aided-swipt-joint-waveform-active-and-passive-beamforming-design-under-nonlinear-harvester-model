clear; clc; setup; config_re_irs;

%% ! R-E region for ideal, adaptive, nonadaptive and no IRS
reIdealIrsSample = cell(nChannels, length(Variable.bandwidth));
reAdaptiveIrsSample = cell(nChannels, length(Variable.bandwidth));
reWitIrsSample = cell(nChannels, length(Variable.bandwidth));
reWptIrsSample = cell(nChannels, length(Variable.bandwidth));
reNoIrsSample = cell(nChannels, length(Variable.bandwidth));

reIdealIrsSolution = cell(nChannels, length(Variable.bandwidth));
reAdaptiveIrsSolution = cell(nChannels, length(Variable.bandwidth));
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
		[directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
		[incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
		[reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
		[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

		% * Upper bound by frequency-selective IRS
		[idealCompositeChannel] = composite_channel_ideal(directChannel, cascadedChannel);
		[reIdealIrsSample{iChannel, iBandwidth}, reIdealIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, idealCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Adaptive IRS and waveform design
		[reAdaptiveIrsSample{iChannel, iBandwidth}, reAdaptiveIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		witCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{1}.compositeChannel;
		wptCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{end}.compositeChannel;

		% * Waveform optimization with nonadaptive WIT-optimized IRS
		[reWitIrsSample{iChannel, iBandwidth}, reWitIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, witCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization with nonadaptive WPT-optimized IRS
		[reWptIrsSample{iChannel, iBandwidth}, reWptIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, wptCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization without IRS
		[reNoIrsSample{iChannel, iBandwidth}, reNoIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
	end
end

% * Average over channel realizations
reIdealIrsInstance = cell(1, length(Variable.bandwidth));
reAdaptiveIrsInstance = cell(1, length(Variable.bandwidth));
reWitIrsInstance = cell(1, length(Variable.bandwidth));
reWptIrsInstance = cell(1, length(Variable.bandwidth));
reNoIrsInstance = cell(1, length(Variable.bandwidth));

for iBandwidth = 1 : length(Variable.bandwidth)
	reIdealIrsInstance{iBandwidth} = mean(cat(3, reIdealIrsSample{:, iBandwidth}), 3);
	reAdaptiveIrsInstance{iBandwidth} = mean(cat(3, reAdaptiveIrsSample{:, iBandwidth}), 3);
	reWitIrsInstance{iBandwidth} = mean(cat(3, reWitIrsSample{:, iBandwidth}), 3);
	reWptIrsInstance{iBandwidth} = mean(cat(3, reWptIrsSample{:, iBandwidth}), 3);
	reNoIrsInstance{iBandwidth} = mean(cat(3, reNoIrsSample{:, iBandwidth}), 3);
end

% * Save batch data
save(sprintf('data/re_irs/re_irs_%d.mat', iBatch));
