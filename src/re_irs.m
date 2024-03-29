clear; clc; setup; config_re_irs;

%% ! R-E region for ideal, adaptive, nonadaptive and no IRS
reIdealIrsSample = cell(nChannels, length(Variable.bandwidth));
reAdaptiveIrsSample = cell(nChannels, length(Variable.bandwidth));
reLinearIrsSample = cell(nChannels, length(Variable.bandwidth));
reWitIrsSample = cell(nChannels, length(Variable.bandwidth));
reWptIrsSample = cell(nChannels, length(Variable.bandwidth));
reRandomIrsSample = cell(nChannels, length(Variable.bandwidth));
reNoIrsSample = cell(nChannels, length(Variable.bandwidth));

reIdealIrsSolution = cell(nChannels, length(Variable.bandwidth));
reAdaptiveIrsSolution = cell(nChannels, length(Variable.bandwidth));
reLinearIrsSolution = cell(nChannels, length(Variable.bandwidth));
reWitIrsSolution = cell(nChannels, length(Variable.bandwidth));
reWptIrsSolution = cell(nChannels, length(Variable.bandwidth));
reRandomIrsSolution = cell(nChannels, length(Variable.bandwidth));
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
		[reIdealIrsSample{iChannel, iBandwidth}, reIdealIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_benchmark(alpha, beta2, beta4, idealCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Adaptive IRS and waveform design
		[reAdaptiveIrsSample{iChannel, iBandwidth}, reAdaptiveIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		witCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{1}.compositeChannel;
		wptCompositeChannel = reAdaptiveIrsSolution{iChannel, iBandwidth}{end}.compositeChannel;

		% * IRS design based on linear EH model, waveform design based on nonlinear model
		[reLinearIrsSample{iChannel, iBandwidth}, reLinearIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_linear_irs(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

		% * Waveform optimization with nonadaptive WIT-optimized IRS
		[reWitIrsSample{iChannel, iBandwidth}, reWitIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_benchmark(alpha, beta2, beta4, witCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization with nonadaptive WPT-optimized IRS
		[reWptIrsSample{iChannel, iBandwidth}, reWptIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_benchmark(alpha, beta2, beta4, wptCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform design with random IRS
		randomIrs = exp(1i * 2 * pi * rand(nReflectors, 1));
		randomCompositeChannel = composite_channel(directChannel, cascadedChannel, randomIrs);
		[reRandomIrsSample{iChannel, iBandwidth}, reRandomIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_benchmark(alpha, beta2, beta4, randomCompositeChannel, txPower, noisePower, nSamples, tolerance);

		% * Waveform optimization without IRS
		[reNoIrsSample{iChannel, iBandwidth}, reNoIrsSolution{iChannel, iBandwidth}] = re_sample_swipt_benchmark(alpha, beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
	end
end

% * Average over channel realizations
reIdealIrsInstance = cell(1, length(Variable.bandwidth));
reAdaptiveIrsInstance = cell(1, length(Variable.bandwidth));
reLinearIrsInstance = cell(1, length(Variable.bandwidth));
reWitIrsInstance = cell(1, length(Variable.bandwidth));
reWptIrsInstance = cell(1, length(Variable.bandwidth));
reRandomIrsInstance = cell(1, length(Variable.bandwidth));
reNoIrsInstance = cell(1, length(Variable.bandwidth));
flag = zeros(1, length(Variable.bandwidth));

for iBandwidth = 1 : length(Variable.bandwidth)
	reIdealIrsInstance{iBandwidth} = mean(cat(3, reIdealIrsSample{:, iBandwidth}), 3);
	reAdaptiveIrsInstance{iBandwidth} = mean(cat(3, reAdaptiveIrsSample{:, iBandwidth}), 3);
	reLinearIrsInstance{iBandwidth} = mean(cat(3, reLinearIrsSample{:, iBandwidth}), 3);
	reWitIrsInstance{iBandwidth} = mean(cat(3, reWitIrsSample{:, iBandwidth}), 3);
	reWptIrsInstance{iBandwidth} = mean(cat(3, reWptIrsSample{:, iBandwidth}), 3);
	reRandomIrsInstance{iBandwidth} = mean(cat(3, reRandomIrsSample{:, iBandwidth}), 3);
	reNoIrsInstance{iBandwidth} = mean(cat(3, reNoIrsSample{:, iBandwidth}), 3);
	flag(iBandwidth) = isempty(reIdealIrsInstance{iBandwidth}) || isempty(reAdaptiveIrsInstance{iBandwidth}) || isempty(reLinearIrsInstance{iBandwidth}) || isempty(reWitIrsInstance{iBandwidth}) || isempty(reWptIrsInstance{iBandwidth}) || isempty(reRandomIrsInstance{iBandwidth}) || isempty(reNoIrsInstance{iBandwidth});
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_irs/re_irs_%d.mat', iBatch));
end
