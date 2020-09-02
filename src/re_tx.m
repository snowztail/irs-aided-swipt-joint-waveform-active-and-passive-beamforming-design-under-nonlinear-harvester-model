clear; clc; setup; config_tx;

%% ! R-E region vs number of transmit antennas
reSample = cell(nChannels, length(Variable.nTxs));
reSolution = cell(nChannels, length(Variable.nTxs));

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
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

        % * Alternating optimization
        [reSample{iChannel, iTx}, reSolution{iChannel, iTx}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reSampleAvg = cell(1, length(Variable.nTxs));
for iTx = 1 : length(Variable.nTxs)
    reSampleAvg{iTx} = mean(cat(3, reSample{:, iTx}), 3);
end

% * Save data
load('data/re_tx.mat');
reSet(iBatch, :) = reSampleAvg;
save('data/re_tx.mat', 'reSet', '-append');
