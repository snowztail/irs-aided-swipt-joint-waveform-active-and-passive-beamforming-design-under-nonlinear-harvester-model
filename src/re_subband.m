clear; clc; setup; config_re_subband;

%% ! R-E region vs number of subbands
reAoSample = cell(nChannels, length(Variable.nSubbands));
reLcSample = cell(nChannels, length(Variable.nSubbands));

reAoSolution = cell(nChannels, length(Variable.nSubbands));
reLcSolution = cell(nChannels, length(Variable.nSubbands));

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
        [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
		[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

        % * Alternating optimization
		[reAoSample{iChannel, iSubband}, reAoSolution{iChannel, iSubband}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		[reLcSample{iChannel, iSubband}, reLcSolution{iChannel, iSubband}] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAoInstance = cell(1, length(Variable.nSubbands));
reLcInstance = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
	reAoInstance{iSubband} = mean(cat(3, reAoSample{:, iSubband}), 3);
	reLcInstance{iSubband} = mean(cat(3, reLcSample{:, iSubband}), 3);
end

% * Save batch data
save(sprintf('data/re_subband/re_subband_%d.mat', iBatch));
