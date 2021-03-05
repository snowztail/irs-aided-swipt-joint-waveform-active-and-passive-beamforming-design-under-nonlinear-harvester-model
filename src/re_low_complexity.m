clear; clc; setup; config_re_low_complexity;

%% ! R-E region for IRS-aided NLoS and LoS channels
reLcSample = cell(nChannels, 1);
reGpSample = cell(nChannels, 1);

reLcSolution = cell(nChannels, 1);
reGpSolution = cell(nChannels, 1);

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
    [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

	for iAlpha = 1 : length(Variable.alpha)
		% * Get alpha for SMF
		alpha = Variable.alpha(iAlpha);

		% * Optimize IRS for given waveform
		[reLcSample, reLcSolution] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
	end
end

% * Average over channel realizations
reNlosInstance = mean(cat(3, reLcSample{:}), 3);
reLosInstance = mean(cat(3, reGpSample{:}), 3);

% * Save batch data
save(sprintf('data/re_los/re_los_%d.mat', iBatch));
