clear; clc; setup; config_scaling_tx;

%% ! WIT/WPT vs number of IRS elements
witSample = zeros(nChannels, length(Variable.nTxs), 2);
wptLinearSample = zeros(nChannels, length(Variable.nTxs), 2);
wptNonlinearSample = zeros(nChannels, length(Variable.nTxs), 2);

witSolution = cell(nChannels, length(Variable.nTxs));
wptLinearSolution = cell(nChannels, length(Variable.nTxs));
wptNonlinearSolution = cell(nChannels, length(Variable.nTxs));

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

		% * WIT
		[witSample(iChannel, iTx, :), witSolution{iChannel, iTx}] = re_sample_wit_wf(directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by linear harvester model
		[wptLinearSample(iChannel, iTx, :), wptLinearSolution{iChannel, iTx}] = re_sample_wpt_ass(beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by nonlinear harvester model
		[wptNonlinearSample(iChannel, iTx, :), wptNonlinearSolution{iChannel, iTx}] = re_sample_wpt_sdr(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);
    end
end

% * Average over channel realizations
rateInstance = mean(witSample(:, :, 1), 1);
currentLinearInstance = mean(wptLinearSample(:, :, 2), 1);
currentNonlinearInstance = mean(wptNonlinearSample(:, :, 2), 1);

% * Save batch data
save(sprintf('data/scaling_tx/scaling_tx_%d.mat', iBatch));
