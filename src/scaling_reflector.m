clear; clc; setup; config_scaling_reflector;

%% ! WIT/WPT vs number of IRS elements
witWfSample = zeros(nChannels, length(Variable.nReflectors), 2);
wptAssSample = zeros(nChannels, length(Variable.nReflectors), 2);
wptSmfSample = zeros(nChannels, length(Variable.nReflectors), 2);
wptSdrSample = zeros(nChannels, length(Variable.nReflectors), 2);

witWfSolution = cell(nChannels, length(Variable.nReflectors));
wptAssSolution = cell(nChannels, length(Variable.nReflectors));
wptSmfSolution = cell(nChannels, length(Variable.nReflectors));
wptSdrSolution = cell(nChannels, length(Variable.nReflectors));

for iChannel = 1 : nChannels
    for iReflector = 1 : length(Variable.nReflectors)
        % * Get number of reflectors and define spatial correlation
        nReflectors = Variable.nReflectors(iReflector);
        corIrs = eye(nReflectors);

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
		[witWfSample(iChannel, iReflector, :), witWfSolution{iChannel, iReflector}] = re_sample_wit_wf(directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by adaptive single sine
		[wptAssSample(iChannel, iReflector, :), wptAssSolution{iChannel, iReflector}] = re_sample_wpt_ass(beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by scaled matched filter
		[wptSmfSample(iChannel, iReflector, :), wptSmfSolution{iChannel, iReflector}] = re_sample_wpt_smf(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by semi-definite relaxation
		[wptSdrSample(iChannel, iReflector, :), wptSdrSolution{iChannel, iReflector}] = re_sample_wpt_sdr(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance);
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
	save(sprintf('data/scaling_reflector/scaling_reflector_%d.mat', iBatch));
end
