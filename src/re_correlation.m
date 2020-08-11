clear; clc; setup; config_correlation;

%% ! R-E region for spatial correlated and uncorrelated channels
% * Uncorrelated channels
[directUncorrelatedChannel] = frequency_response(directUncorrelatedTapGain, directTapDelay, directDistance, subbandFrequency, fadingMode);
[incidentUncorrelatedChannel] = frequency_response(incidentUncorrelatedTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
[reflectiveUncorrelatedChannel] = frequency_response(reflectiveUncorrelatedTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);
[reSample{1}, reSolution{1}] = re_sample(beta2, beta4, directUncorrelatedChannel, incidentUncorrelatedChannel, reflectiveUncorrelatedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

% * Correlated channels
[directCorrelatedChannel] = frequency_response(directCorrelatedTapGain, directTapDelay, directDistance, subbandFrequency, fadingMode);
[incidentCorrelatedChannel] = frequency_response(incidentCorrelatedTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
[reflectiveCorrelatedChannel] = frequency_response(reflectiveCorrelatedTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);
[reSample{2}, reSolution{2}] = re_sample(beta2, beta4, directCorrelatedChannel, incidentCorrelatedChannel, reflectiveCorrelatedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

save('data/re_correlation.mat');

%% * R-E plots
figure('name', 'R-E region for spatial correlated and uncorrelated channels');
plot(reSample{1}(1, :) / nSubbands, 1e6 * reSample{1}(2, :));
hold on;
plot(reSample{2}(1, :) / nSubbands, 1e6 * reSample{2}(2, :));
hold off;
grid minor;
legend('Uncorrelated', 'Correlated');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_correlation.fig');
