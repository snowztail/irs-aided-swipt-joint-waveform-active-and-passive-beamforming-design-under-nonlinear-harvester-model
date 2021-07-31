clear; clc; setup; config_re_lc;

%% ! R-E region vs SMF ratio
reAoSample = cell(nChannels, length(Variable.nSubbands));
reLcSample = cell(nChannels, length(Variable.nSubbands), length(Variable.alpha));

reAoSolution = cell(nChannels, length(Variable.nSubbands));
reLcSolution = cell(nChannels, length(Variable.nSubbands), length(Variable.alpha));

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

        % * R-E region by AO
		[reAoSample{iChannel, iSubband}, reAoSolution{iChannel, iSubband}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

		% * R-E region by low-complexity design
		for iAlpha = 1 : length(Variable.alpha)
			[reLcSample{iChannel, iSubband, iAlpha}, reLcSolution{iChannel, iSubband, iAlpha}] = re_sample_swipt_low_complexity(Variable.alpha(iAlpha), beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		end
    end
end

% * Average over channel realizations
reAoInstance = cell(length(Variable.nSubbands), 1);
reLcInstance = cell(length(Variable.nSubbands), length(Variable.alpha));
flag = zeros(length(Variable.nSubbands), length(Variable.alpha));

for iSubband = 1 : length(Variable.nSubbands)
	reAoInstance{iSubband} = mean(cat(3, reAoSample{:, iSubband}), 3);
	for iAlpha = 1 : length(Variable.alpha)
		reLcInstance{iSubband, iAlpha} = mean(cat(4, reLcSample{:, iSubband, iAlpha}), 4);
		flag(iSubband, iAlpha) = isempty(reAoInstance{iSubband}) || isempty(reLcInstance{iSubband, iAlpha});
	end
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_lc/re_lc_%d.mat', iBatch));
end
