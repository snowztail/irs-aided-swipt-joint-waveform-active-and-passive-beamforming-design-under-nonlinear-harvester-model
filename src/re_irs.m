clear; clc; setup; config_irs;

%% ! R-E region for fixed and adaptive IRS
% * Generate channels
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, subbandFrequency, fadingMode);
[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);

% * Adaptive IRS and waveform design
[reSample{1}, reSolution{1}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
compositeChannelWit = reSolution{1}{1}.compositeChannel;
compositeChannelWpt = reSolution{1}{end}.compositeChannel;

% * Waveform optimization with fixed WIT-optimized IRS
[reSample{2}, reSolution{2}] = re_sample_reference(beta2, beta4, compositeChannelWit, txPower, noisePower, nSamples, tolerance);

% * Waveform optimization with fixed WPT-optimized IRS
[reSample{3}, reSolution{3}] = re_sample_reference(beta2, beta4, compositeChannelWpt, txPower, noisePower, nSamples, tolerance);

% * Waveform optimization without IRS
[reSample{4}, reSolution{4}] = re_sample_reference(beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);

save('data/re_irs.mat');

%% * R-E plots
figure('name', 'R-E region for adaptive, fixed and no IRS');
plot(reSample{1}(1, :) / nSubbands, 1e6 * reSample{1}(2, :));
hold on;
plot(reSample{2}(1, :) / nSubbands, 1e6 * reSample{2}(2, :));
hold on;
plot(reSample{3}(1, :) / nSubbands, 1e6 * reSample{3}(2, :));
hold on;
plot(reSample{4}(1, :) / nSubbands, 1e6 * reSample{4}(2, :));
hold off;
grid minor;
legend('Adaptive IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_irs.fig');
