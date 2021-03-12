clear; clc; setup; config_re_csi;

%% ! R-E region vs CSIT error
reNoIrsSample = cell(nChannels, 1);
reErrorSample = cell(nChannels, length(Variable.cascadedErrorVariance));

reNoIrsSolution = cell(nChannels, 1);
reErrorSolution = cell(nChannels, length(Variable.cascadedErrorVariance));

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

	% * R-E region without IRS
	[reNoIrsSample{iChannel}, reNoIrsSolution{iChannel}] = re_sample_swipt_gp_benchmark(alpha, beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);

    for iError = 1 : length(Variable.cascadedErrorVariance)
		% * Update error variance of the cascaded channel
        cascadedErrorVariance = Variable.cascadedErrorVariance(iError);

		% * Emulate cascaded channel estimation
		[imperfectCascadedChannel] = imperfect_csi(cascadedChannel, cascadedErrorVariance);

        % * R-E region based on imperfect cascaded CSIT
        [reErrorSample{iChannel, iError}, reErrorSolution{iChannel, iError}] = re_sample_swipt_imperfect_csi(alpha, beta2, beta4, directChannel, cascadedChannel, imperfectCascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reNoIrsInstance = mean(cat(3, reNoIrsSample{:}), 3);
reErrorInstance = cell(1, length(Variable.cascadedErrorVariance));
for iError = 1 : length(Variable.cascadedErrorVariance)
    reErrorInstance{iError} = mean(cat(3, reErrorSample{:, iError}), 3);
end

% * Save batch data
save(sprintf('data/re_csi/re_csi_%d.mat', iBatch));
