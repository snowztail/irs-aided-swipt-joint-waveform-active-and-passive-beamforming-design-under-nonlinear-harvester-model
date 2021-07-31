clear; clc; setup; config_scaling_tx;

%% ! WIT/WPT vs number of IRS elements
witWfSample = zeros(nChannels, length(Variable.nTxs), 2);
wptAssSample = zeros(nChannels, length(Variable.nTxs), 2);
wptSmfSample = zeros(nChannels, length(Variable.nTxs), 2);
wptSdrSample = zeros(nChannels, length(Variable.nTxs), 2);

witWfSolution = cell(nChannels, length(Variable.nTxs));
wptAssSolution = cell(nChannels, length(Variable.nTxs));
wptSmfSolution = cell(nChannels, length(Variable.nTxs));
wptSdrSolution = cell(nChannels, length(Variable.nTxs));

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

		% * WIT by water filling
		[witWfSample(iChannel, iTx, :), witWfSolution{iChannel, iTx}] = re_sample_wit_wf(directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by adaptive single sine
		[wptAssSample(iChannel, iTx, :), wptAssSolution{iChannel, iTx}] = re_sample_wpt_ass(beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by scaled matched filter
		[wptSmfSample(iChannel, iTx, :), wptSmfSolution{iChannel, iTx}] = re_sample_wpt_smf(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by semi-definite relaxation
		[wptSdrSample(iChannel, iTx, :), wptSdrSolution{iChannel, iTx}] = re_sample_wpt_sdr(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);
    end
end

% * Average over channel realizations
rateWf = mean(witWfSample(:, :, 1), 1);
currentAss = mean(wptAssSample(:, :, 2), 1);
currentSmf = mean(wptSmfSample(:, :, 2), 1);
currentSdr = mean(wptSdrSample(:, :, 2), 1);
flag = isempty(rateWf) || isempty(currentAss) || isempty(currentSmf) || isempty(currentSdr);

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/scaling_tx/scaling_tx_%d.mat', iBatch));
end
