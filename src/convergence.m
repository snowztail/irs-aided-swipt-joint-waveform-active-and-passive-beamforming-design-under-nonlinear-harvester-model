clear; clc; setup; config_convergence;

%% ! Convergence test for all algorithms

% * Generate tap gains and delays
[directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
[incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

% * Construct channels
[directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
[incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
[reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

[aoSample, aoSolution] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
[lcSample, lcSolution] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

flag = isempty(aoSample) || isempty(lcSample);

% * Save batch data
if ~sum(flag(:))
	save('data/convergence.mat');
end
