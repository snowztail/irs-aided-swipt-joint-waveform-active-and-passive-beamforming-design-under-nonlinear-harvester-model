clear; clc; setup; config_scaling_reflector;

%% ! WIT/WPT vs number of IRS elements
witSample = zeros(nChannels, length(Variable.nReflectors), 2);
wptLinearSample = zeros(nChannels, length(Variable.nReflectors), 2);
wptNonlinearSample = zeros(nChannels, length(Variable.nReflectors), 2);

witSolution = cell(nChannels, length(Variable.nReflectors));
wptLinearSolution = cell(nChannels, length(Variable.nReflectors));
wptNonlinearSolution = cell(nChannels, length(Variable.nReflectors));

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

		% * WIT
		[witSample(iChannel, iReflector, :), witSolution{iChannel, iReflector}] = re_sample_wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by linear harvester model
		[wptLinearSample(iChannel, iReflector, :), wptLinearSolution{iChannel, iReflector}] = re_sample_wpt_linear(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);

		% * WPT by nonlinear harvester model
		[wptNonlinearSample(iChannel, iReflector, :), wptNonlinearSolution{iChannel, iReflector}] = re_sample_wpt(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
    end
end

% * Average over channel realizations
rateInstance = mean(witSample(:, :, 1), 1);
currentLinearInstance = mean(wptLinearSample(:, :, 2), 1);
currentNonlinearInstance = mean(wptNonlinearSample(:, :, 2), 1);

% * Save batch data
save(sprintf('data/scaling_reflector/scaling_reflector_%d.mat', iBatch));
