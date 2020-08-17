clear; clc; setup; config_distance;

%% ! R-E region vs AP-IRS distance
reSample = cell(nChannels, nCases);
reSolution = cell(nChannels, nCases);

parfor iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct direct channel
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);

    for iDistance = 1 : nCases
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
reSampleAvg = cell(1, nCases);
for iDistance = 1 : nCases
    reSampleAvg{iDistance} = mean(cat(3, reSample{:, iDistance}), 3);
end
save('data/re_distance.mat');

% %% * R-E plots
% figure('name', 'R-E region vs AP-IRS horizontal distance');
% legendString = cell(1, nCases);
% for iDistance = 1 : nCases
%     plot(reSampleAvg{iDistance}(1, :) / nSubbands, 1e6 * reSampleAvg{iDistance}(2, :));
%     legendString{iDistance} = sprintf('d_H = %d', Variable.horizontalDistance(iDistance));
%     hold on;
% end
% hold off;
% grid minor;
% legend(legendString);
% xlabel('Per-subband rate [bps/Hz]');
% ylabel('Average output DC current [\muA]');
% ylim([0 inf]);
% savefig('plots/re_distance.fig');
