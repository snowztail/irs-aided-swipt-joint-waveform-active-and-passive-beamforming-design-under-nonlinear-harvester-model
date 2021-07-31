clear; clc; setup; config_re_noise;

%% ! R-E region for different large-scale SNR
reAoSample = cell(nChannels, length(Variable.noisePower));
reLcSample = cell(nChannels, length(Variable.noisePower));

reAoSolution = cell(nChannels, length(Variable.noisePower));
reLcSolution = cell(nChannels, length(Variable.noisePower));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
    [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
	[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

    for iNoise = 1 : length(Variable.noisePower)
        % * Get average noise power
        noisePower = Variable.noisePower(iNoise);

        % * Alternating optimization
        [reAoSample{iChannel, iNoise}, reAoSolution{iChannel, iNoise}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
        [reLcSample{iChannel, iNoise}, reLcSolution{iChannel, iNoise}] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAoInstance = cell(1, length(Variable.noisePower));
reLcInstance = cell(1, length(Variable.noisePower));
flag = zeros(1, length(Variable.noisePower));

for iNoise = 1 : length(Variable.noisePower)
    reAoInstance{iNoise} = mean(cat(3, reAoSample{:, iNoise}), 3);
    reLcInstance{iNoise} = mean(cat(3, reLcSample{:, iNoise}), 3);
	flag(iNoise) = isempty(reAoInstance{iNoise}) || isempty(reLcInstance{iNoise});
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_noise/re_noise_%d.mat', iBatch));
end
