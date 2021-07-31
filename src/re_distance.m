clear; clc; setup; config_re_distance;

%% ! R-E region vs AP-IRS distance
reAoSample = cell(nChannels, length(Variable.horizontalDistance));
reLcSample = cell(nChannels, length(Variable.horizontalDistance));

reAoSolution = cell(nChannels, length(Variable.horizontalDistance));
reLcSolution = cell(nChannels, length(Variable.horizontalDistance));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct direct channel
    [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);

    for iDistance = 1 : length(Variable.horizontalDistance)
        % * Get distances
        horizontalDistance = Variable.horizontalDistance(iDistance);
        [incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance);

        % * Construct extra channels
        [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
		cascadedChannel = cascaded_channel(incidentChannel, reflectiveChannel);

        % * Alternating optimization
        [reAoSample{iChannel, iDistance}, reAoSolution{iChannel, iDistance}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		[reLcSample{iChannel, iDistance}, reLcSolution{iChannel, iDistance}] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAoInstance = cell(1, length(Variable.horizontalDistance));
reLcInstance = cell(1, length(Variable.horizontalDistance));
flag = zeros(1, length(Variable.horizontalDistance));

for iDistance = 1 : length(Variable.horizontalDistance)
    reAoInstance{iDistance} = mean(cat(3, reAoSample{:, iDistance}), 3);
    reLcInstance{iDistance} = mean(cat(3, reLcSample{:, iDistance}), 3);
	flag(iDistance) = isempty(reAoInstance{iDistance}) || isempty(reLcInstance{iDistance});
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_distance/re_distance_%d.mat', iBatch));
end
