clear; clc; setup; config_irs; load('data/tap_los.mat');

%% ! R-E region for fixed and adaptive IRS
reSample = cell(2, 1);
reSolution = cell(2, 1);

% * Get tap data
directTapGain = directTapGain(:, 1 : nTxs);
incidentTapGain = incidentTapGain(:, 1 : nTxs, :);

% * Generate channels
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

% * Waveform optimization (with fixed WIT-optimized IRS)
[reSample{1}, reSolution{1}] = re_sample_fixed_irs(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

% * Alternating optimization
[reSample{2}, reSolution{2}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

save('data/re_irs.mat');

%% * R-E plots
figure('name', 'R-E region for fixed and adaptive IRS');
plot(reSample{1}(1, :) / nSubbands, 1e6 * reSample{1}(2, :));
hold on;
plot(reSample{2}(1, :) / nSubbands, 1e6 * reSample{2}(2, :));
hold off;
grid minor;
legend('Fixed IRS', 'Adaptive IRS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_irs.fig');
