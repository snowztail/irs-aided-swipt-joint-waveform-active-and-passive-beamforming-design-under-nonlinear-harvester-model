clear; clc; setup; config_ff_gp_sdr; load('data/tap.mat');

%% ! FF-IRS: R-E region by GP and SDR
% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

% * GP and SDR
ff_gp;
ff_sdr;
save('data/re_ff_gp_sdr.mat');

%% * R-E plots
figure('name', 'FF-IRS: R-E region by GP and SDR')
plot(ffGpSample(1, :) / nSubbands, 1e6 * ffGpSample(2, :), 'k--');
hold on;
plot(ffSdrSample(1, :) / nSubbands, 1e6 * ffSdrSample(2, :), 'r-');
hold off;
grid minor;
legend('GP', 'SDR');
xlabel('Rate per subband [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_ff_gp_sdr.fig');
