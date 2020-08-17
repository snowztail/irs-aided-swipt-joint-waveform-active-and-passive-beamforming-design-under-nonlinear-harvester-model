clear; clc; setup; config_irs;

%% ! R-E region for fixed and adaptive IRS
reSample = cell(nChannels, nCases);
reSolution = cell(nChannels, nCases);

parfor iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
    [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    % * Adaptive IRS and waveform design
    [reSample{iChannel, 1}, reSolution{iChannel, 1}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    compositeChannelWit = reSolution{iChannel, 1}{1}.compositeChannel;
    compositeChannelWpt = reSolution{iChannel, 1}{end}.compositeChannel;

    % * Waveform optimization with fixed WIT-optimized IRS
    [reSample{iChannel, 2}, reSolution{iChannel, 2}] = re_sample_reference(beta2, beta4, compositeChannelWit, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization with fixed WPT-optimized IRS
    [reSample{iChannel, 3}, reSolution{iChannel, 3}] = re_sample_reference(beta2, beta4, compositeChannelWpt, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization without IRS
    [reSample{iChannel, 4}, reSolution{iChannel, 4}] = re_sample_reference(beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
end

% * Average over channel realizations
reSampleAvg = cell(1, nCases);
for iCase = 1 : nCases
    reSampleAvg{iCase} = mean(cat(3, reSample{:, iCase}), 3);
end
save('data/re_irs.mat');

% %% * R-E plots
% figure('name', 'R-E region for adaptive, fixed and no IRS');
% for iCase = 1 : nCases
%     plot(reSampleAvg{iCase}(1, :) / nSubbands, 1e6 * reSampleAvg{iCase}(2, :));
%     hold on;
% end
% hold off;
% grid minor;
% legend('Adaptive IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
% xlabel('Per-subband rate [bps/Hz]');
% ylabel('Average output DC current [\muA]');
% ylim([0 inf]);
% savefig('plots/re_irs.fig');
