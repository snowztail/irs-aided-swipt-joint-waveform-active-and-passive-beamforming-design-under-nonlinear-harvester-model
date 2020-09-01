clear; clc; setup; config_irs;

%% ! R-E region for fixed and adaptive IRS
reAdaptiveSample = cell(nChannels, 1);
reWitSample = cell(nChannels, 1);
reWptSample = cell(nChannels, 1);
reNoIrsSample = cell(nChannels, 1);
reAdaptiveSolution = cell(nChannels, 1);
reWitSolution = cell(nChannels, 1);
reWptSolution = cell(nChannels, 1);
reNoIrsSolution = cell(nChannels, 1);

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
    [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    % * Adaptive IRS and waveform design
    [reAdaptiveSample{iChannel}, reAdaptiveSolution{iChannel}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    compositeChannelWit = reAdaptiveSolution{iChannel}{1}.compositeChannel;
    compositeChannelWpt = reAdaptiveSolution{iChannel}{end}.compositeChannel;

    % * Waveform optimization with fixed WIT-optimized IRS
    [reWitSample{iChannel}, reWitSolution{iChannel}] = re_sample_reference(beta2, beta4, compositeChannelWit, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization with fixed WPT-optimized IRS
    [reWptSample{iChannel}, reWptSolution{iChannel}] = re_sample_reference(beta2, beta4, compositeChannelWpt, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization without IRS
    [reNoIrsSample{iChannel}, reNoIrsSolution{iChannel}] = re_sample_reference(beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
end

% * Average over channel realizations
reAdaptiveSampleAvg = mean(cat(3, reAdaptiveSample{:}), 3);
reWitSampleAvg = mean(cat(3, reWitSample{:}), 3);
reWptSampleAvg = mean(cat(3, reWptSample{:}), 3);
reNoIrsSampleAvg = mean(cat(3, reNoIrsSample{:}), 3);
reSampleAvg = [reAdaptiveSampleAvg; reWitSampleAvg; reWptSampleAvg; reNoIrsSampleAvg];

% * Save data
load('data/re_irs.mat');
reSet(:, pbsIndex) = reSampleAvg;
save('data/re_irs.mat', 'reSet', '-append');

% %% * R-E plots
% figure('name', 'R-E region for adaptive, fixed and no IRS');
% plot(reAdaptiveSampleAvg(1, :) / nSubbands, 1e6 * reAdaptiveSampleAvg(2, :));
% hold on;
% plot(reWitSampleAvg(1, :) / nSubbands, 1e6 * reWitSampleAvg(2, :));
% hold on;
% plot(reWptSampleAvg(1, :) / nSubbands, 1e6 * reWptSampleAvg(2, :));
% hold on;
% plot(reNoIrsSampleAvg(1, :) / nSubbands, 1e6 * reNoIrsSampleAvg(2, :));
% hold off;
% grid minor;
% legend('Adaptive IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
% xlabel('Per-subband rate [bps/Hz]');
% ylabel('Average output DC current [\muA]');
% ylim([0 inf]);
% savefig('plots/re_irs.fig');
