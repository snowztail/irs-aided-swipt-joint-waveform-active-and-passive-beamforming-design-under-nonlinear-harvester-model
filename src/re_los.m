clear; clc; setup; config_los;

%% ! R-E region for IRS-aided NLOS and LOS channels
reNlosSample = cell(nChannels, 1);
reLosSample = cell(nChannels, 1);
reNlosSolution = cell(nChannels, 1);
reLosSolution = cell(nChannels, 1);

parfor iChannel = 1 : nChannels
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

    % * Optimization based on NLOS channels
    [reNlosSample{iChannel}, reNlosSolution{iChannel}] = re_sample(beta2, beta4, directNlosChannel, incidentNlosChannel, reflectiveNlosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

    % * Optimization based on NLOS direct channel and LOS incident and reflective channels
    [reLosSample{iChannel}, reLosSolution{iChannel}] = re_sample(beta2, beta4, directNlosChannel, incidentLosChannel, reflectiveLosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end

% * Average over channel realizations
reNlosSampleAvg = mean(cat(3, reNlosSample{:}), 3);
reLosSampleAvg = mean(cat(3, reLosSample{:}), 3);
save('data/re_los.mat');

% %% * R-E plots
% figure('name', 'R-E region for IRS-aided NLOS and LOS channels');
% plot(reNlosSampleAvg(1, :) / nSubbands, 1e6 * reNlosSampleAvg(2, :));
% hold on;
% plot(reLosSampleAvg(1, :) / nSubbands, 1e6 * reLosSampleAvg(2, :));
% hold off;
% grid minor;
% legend('NLOS', 'LOS');
% xlabel('Per-subband rate [bps/Hz]');
% ylabel('Average output DC current [\muA]');
% ylim([0 inf]);
% savefig('plots/re_los.fig');
