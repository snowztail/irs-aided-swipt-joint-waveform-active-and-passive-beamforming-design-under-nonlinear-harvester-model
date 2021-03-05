clear; clc; setup; config_re_subband;

%% ! R-E region vs number of subbands
reSample = cell(nChannels, length(Variable.nSubbands));
reSolution = cell(nChannels, length(Variable.nSubbands));

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

        % * Alternating optimization
		[reSample{iChannel, iSubband}, reSolution{iChannel, iSubband}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reInstance = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
	reInstance{iSubband} = mean(cat(3, reSample{:, iSubband}), 3);
end

% * Save batch data
save(sprintf('data/re_subband/re_subband_%d.mat', iBatch));
