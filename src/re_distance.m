clear; clc; setup; config_distance;

%% ! R-E region vs AP-IRS distance
reSample = cell(nChannels, length(Variable.horizontalDistance));
reSolution = cell(nChannels, length(Variable.horizontalDistance));

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct direct channel
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);

    for iDistance = 1 : length(Variable.horizontalDistance)
        % * Get distances
        horizontalDistance = Variable.horizontalDistance(iDistance);
        [incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance);

        % * Construct extra channels
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

        % * Alternating optimization
        [reSample{iChannel, iDistance}, reSolution{iChannel, iDistance}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reSampleAvg = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    reSampleAvg{iDistance} = mean(cat(3, reSample{:, iDistance}), 3);
end

% * Save data
load('data/re_distance.mat');
reSet(iBatch, :) = reSampleAvg;
save('data/re_distance.mat', 'reSet', '-append');
