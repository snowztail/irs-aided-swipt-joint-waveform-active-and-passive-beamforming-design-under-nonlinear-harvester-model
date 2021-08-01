clear; clc; setup; config_re_waveform;

%% ! WIT/WPT waveform with and without IRS vs number of subbands
reIrsSample = cell(nChannels, length(Variable.nSubbands));
reNoIrsSample = cell(nChannels, length(Variable.nSubbands));

reIrsSolution = cell(nChannels, length(Variable.nSubbands));
reNoIrsSolution = cell(nChannels, length(Variable.nSubbands));

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

		% * IRS
		[reIrsSample{iChannel, iSubband}, reIrsSolution{iChannel, iSubband}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

		% * No IRS
		[reNoIrsSample{iChannel, iSubband}, reNoIrsSolution{iChannel, iSubband}] = re_sample_swipt_benchmark(alpha, beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
    end
end

% * Average over channel realizations
reIrsInstance = cell(1, length(Variable.nSubbands));
reNoIrsInstance = cell(1, length(Variable.nSubbands));
flag = zeros(1, length(Variable.nSubbands));

for iSubband = 1 : length(Variable.nSubbands)
	reIrsInstance{iSubband} = mean(cat(3, reIrsSample{:, iSubband}), 3);
	reNoIrsInstance{iSubband} = mean(cat(3, reNoIrsSample{:, iSubband}), 3);
	flag(iSubband) = isempty(reIrsInstance{iSubband}) || isempty(reNoIrsInstance{iSubband});
end

% * Save batch data
if ~sum(flag(:))
	save('data/re_waveform.mat');
end
