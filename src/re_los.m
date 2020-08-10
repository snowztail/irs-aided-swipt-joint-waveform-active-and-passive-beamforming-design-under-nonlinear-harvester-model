clear; clc; setup; config_los;

%% ! R-E region for IRS-aided NLOS and LOS channels
% * NLOS channels
[directNlosChannel] = frequency_response(directNlosTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
[incidentNlosChannel] = frequency_response(incidentNlosTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveNlosChannel] = frequency_response(reflectiveNlosTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');
[reSample{1}, reSolution{1}] = re_sample(beta2, beta4, directNlosChannel, incidentNlosChannel, reflectiveNlosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

% * (IRS-aided) LOS channels
[incidentLosChannel] = frequency_response(incidentLosTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveLosChannel] = frequency_response(reflectiveLosTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');
[reSample{2}, reSolution{2}] = re_sample(beta2, beta4, directNlosChannel, incidentLosChannel, reflectiveLosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

save('data/re_los.mat');

%% * R-E plots
figure('name', 'R-E region for IRS-aided NLOS and LOS channels');
plot(reSample{1}(1, :) / nSubbands, 1e6 * reSample{1}(2, :));
hold on;
plot(reSample{2}(1, :) / nSubbands, 1e6 * reSample{2}(2, :));
hold off;
grid minor;
legend('NLOS', 'LOS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_los.fig');
