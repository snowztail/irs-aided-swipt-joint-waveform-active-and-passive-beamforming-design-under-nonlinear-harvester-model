clear; clc; setup; config_los;

%% ! R-E region for IRS-aided NLoS and LoS channels
reNlosSample = cell(nChannels, 1);
reLosSample = cell(nChannels, 1);
reNlosSolution = cell(nChannels, 1);
reLosSolution = cell(nChannels, 1);

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directNlosTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentNlosTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveNlosTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    [incidentLosTapGain, ~] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveLosTapGain, ~] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directNlosChannel] = frequency_response(directNlosTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentNlosChannel] = frequency_response(incidentNlosTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
    [reflectiveNlosChannel] = frequency_response(reflectiveNlosTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    [incidentLosChannel] = frequency_response(incidentLosTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
    [reflectiveLosChannel] = frequency_response(reflectiveLosTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    % * Optimization based on NLoS channels
    [reNlosSample{iChannel}, reNlosSolution{iChannel}] = re_sample(beta2, beta4, directNlosChannel, incidentNlosChannel, reflectiveNlosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

    % * Optimization based on NLoS direct channel and LoS incident and reflective channels
    [reLosSample{iChannel}, reLosSolution{iChannel}] = re_sample(beta2, beta4, directNlosChannel, incidentLosChannel, reflectiveLosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end

% * Average over channel realizations
reInstance{1} = mean(cat(3, reNlosSample{:}), 3);
reInstance{2} = mean(cat(3, reLosSample{:}), 3);

% * Save batch data
save(sprintf('data/re_los_%d.mat', iBatch));
