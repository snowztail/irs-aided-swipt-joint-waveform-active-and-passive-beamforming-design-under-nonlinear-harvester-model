clear; clc; setup; config_snr;

%% ! R-E region for different large-scale SNR
reSample = cell(nChannels, length(Variable.snr));
reSolution = cell(nChannels, length(Variable.snr));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
    [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    % * Calculate sum pathloss
    [sumPathloss] = sum_pathloss(directDistance, incidentDistance, reflectiveDistance);

    for iSnr = 1 : length(Variable.snr)
        % * Calculate noise power based on SNR
        noisePower = txPower * sumPathloss * rxGain / Variable.snr(iSnr);

        % * Alternating optimization
        [reSample{iChannel, iSnr}, reSolution{iChannel, iSnr}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reInstance = cell(1, length(Variable.snr));
for iSnr = 1 : length(Variable.snr)
    reInstance{iSnr} = mean(cat(3, reSample{:, iSnr}), 3);
end

% * Save batch data
save(sprintf('data/re_snr_%d.mat', iBatch), 'reInstance');
