clear; clc; setup; config_re_tx;

%% ! R-E region vs number of transmit antennas
reAoSample = cell(nChannels, length(Variable.nTxs));
reLcSample = cell(nChannels, length(Variable.nTxs));

reAoSolution = cell(nChannels, length(Variable.nTxs));
reLcSolution = cell(nChannels, length(Variable.nTxs));

for iChannel = 1 : nChannels
    for iTx = 1 : length(Variable.nTxs)
        % * Get number of transmit antennas and define spatial correlation
        nTxs = Variable.nTxs(iTx);
        corTx = eye(nTxs);

        % * Generate tap gains and delays
        [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
        [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
        [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

        % * Construct channels
        [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
		[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

        % * Alternating optimization
        [reAoSample{iChannel, iTx}, reAoSolution{iChannel, iTx}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
        [reLcSample{iChannel, iTx}, reLcSolution{iChannel, iTx}] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAoInstance = cell(1, length(Variable.nTxs));
reLcInstance = cell(1, length(Variable.nTxs));
for iTx = 1 : length(Variable.nTxs)
    reAoInstance{iTx} = mean(cat(3, reAoSample{:, iTx}), 3);
    reLcInstance{iTx} = mean(cat(3, reLcSample{:, iTx}), 3);
end

% * Save batch data
save(sprintf('data/re_tx/re_tx_%d.mat', iBatch));
