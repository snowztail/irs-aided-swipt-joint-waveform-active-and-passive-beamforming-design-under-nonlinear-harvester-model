clear; clc; setup; config_re_csi;

%% ! R-E region vs CSIT error
reRandomSample = cell(nChannels, length(Variable.nReflectors));
reErrorSample = cell(nChannels, length(Variable.nReflectors), length(Variable.cascadedErrorVariance));

reRandomSolution = cell(nChannels, length(Variable.nReflectors));
reErrorSolution = cell(nChannels, length(Variable.nReflectors), length(Variable.cascadedErrorVariance));

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

		% * R-E region with random IRS
		randomIrs = exp(1i * 2 * pi * rand(nReflectors, 1));
		randomCompositeChannel = composite_channel(directChannel, cascadedChannel, randomIrs);
		[reRandomSample{iChannel, iReflector}, reRandomSolution{iChannel, iReflector}] = re_sample_swipt_benchmark(alpha, beta2, beta4, randomCompositeChannel, txPower, noisePower, nSamples, tolerance);

		for iError = 1 : length(Variable.cascadedErrorVariance)
			% * Update error variance of the cascaded channel
			cascadedErrorVariance = Variable.cascadedErrorVariance(iError);

			% * Emulate cascaded channel estimation
			[imperfectCascadedChannel] = imperfect_csi(cascadedChannel, cascadedErrorVariance);

			% * R-E region based on imperfect cascaded CSIT
			[reErrorSample{iChannel, iReflector, iError}, reErrorSolution{iChannel, iReflector, iError}] = re_sample_swipt_imperfect_csi(alpha, beta2, beta4, directChannel, cascadedChannel, imperfectCascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		end
	end
end

% * Average over channel realizations
reRandomInstance = cell(length(Variable.nReflectors), 1);
reErrorInstance = cell(length(Variable.nReflectors), length(Variable.cascadedErrorVariance));
for iReflector = 1 : length(Variable.nReflectors)
	reRandomInstance{iReflector} = mean(cat(3, reRandomSample{:, iReflector}), 3);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		reErrorInstance{iReflector, iError} = mean(cat(4, reErrorSample{:, iReflector, iError}), 4);
	end
end

% * Save batch data
save(sprintf('data/re_csi/re_csi_%d.mat', iBatch));
