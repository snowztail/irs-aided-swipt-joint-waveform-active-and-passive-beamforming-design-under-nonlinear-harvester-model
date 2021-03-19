clear; clc; setup; config_re_csi;

%% ! R-E region vs CSIT error
reRandomSample = cell(nChannels, length(Variable.nSubbands));
reErrorSample = cell(nChannels, length(Variable.nSubbands), length(Variable.cascadedErrorVariance));

reRandomSolution = cell(nChannels, length(Variable.nSubbands));
reErrorSolution = cell(nChannels, length(Variable.nSubbands), length(Variable.cascadedErrorVariance));

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

		% * R-E region with random IRS
		randomIrs = exp(1i * 2 * pi * rand(nReflectors, 1));
		randomCompositeChannel = composite_channel(directChannel, cascadedChannel, randomIrs);
		[reRandomSample{iChannel, iSubband}, reRandomSolution{iChannel, iSubband}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, randomCompositeChannel, txPower, noisePower, nSamples, tolerance);

		for iError = 1 : length(Variable.cascadedErrorVariance)
			% * Update error variance of the cascaded channel
			cascadedErrorVariance = Variable.cascadedErrorVariance(iError);

			% * Emulate cascaded channel estimation
			[imperfectCascadedChannel] = imperfect_csi(cascadedChannel, cascadedErrorVariance);

			% * R-E region based on imperfect cascaded CSIT
			[reErrorSample{iChannel, iSubband, iError}, reErrorSolution{iChannel, iSubband, iError}] = re_sample_swipt_imperfect_csi(alpha, beta2, beta4, directChannel, cascadedChannel, imperfectCascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
		end
	end
end

% * Average over channel realizations
reRandomInstance = cell(length(Variable.nSubbands), 1);
reErrorInstance = cell(length(Variable.nSubbands), length(Variable.cascadedErrorVariance));
for iSubband = 1 : length(Variable.nSubbands)
	reRandomInstance{iSubband} = mean(cat(3, reRandomSample{:, iSubband}), 3);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		reErrorInstance{iSubband, iError} = mean(cat(4, reErrorSample{:, iSubband, iError}), 4);
	end
end

% * Save batch data
save(sprintf('data/re_csi/re_csi_%d.mat', iBatch));
