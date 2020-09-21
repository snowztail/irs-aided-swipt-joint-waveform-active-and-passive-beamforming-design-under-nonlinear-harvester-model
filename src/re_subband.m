clear; clc; setup; config_subband;

%% ! R-E region vs number of subbands
reAdaptiveIrsSample = cell(nChannels, length(Variable.nSubbands));
reFsIrsSample = cell(nChannels, length(Variable.nSubbands));

reAdaptiveIrsSolution = cell(nChannels, length(Variable.nSubbands));
reFsIrsSolution = cell(nChannels, length(Variable.nSubbands));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    for iSubband = 1 : length(Variable.nSubbands)
        % * Get number of subbands and subband frequency
        nSubbands = Variable.nSubbands(iSubband);
        [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);

        % * Construct channels
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

        % * Adaptive IRS and waveform design
		[reAdaptiveIrsSample{iChannel, iSubband}, reAdaptiveIrsSolution{iChannel, iSubband}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

		% * Upper bound by frequency-selective IRS
		[fsIrsCompositeChannel] = fs_irs_composite_channel(directChannel, incidentChannel, reflectiveChannel);
		[reFsIrsSample{iChannel, iSubband}, reFsIrsSolution{iChannel, iSubband}] = re_sample_reference(beta2, beta4, fsIrsCompositeChannel, txPower, noisePower, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAdaptiveIrsInstance = cell(1, length(Variable.nSubbands));
reFsIrsInstance = cell(1, length(Variable.nSubbands));

for iSubband = 1 : length(Variable.nSubbands)
	reAdaptiveIrsInstance{iSubband} = mean(cat(3, reAdaptiveIrsSample{:, iSubband}), 3);
	reFsIrsInstance{iSubband} = mean(cat(3, reFsIrsSample{:, iSubband}), 3);
end

% * Save batch data
save(sprintf('data/re_subband/re_subband_%d.mat', iBatch));
