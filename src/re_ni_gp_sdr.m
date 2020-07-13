clear; clc; setup; config_ni_gp_sdr;  load('data/tap.mat');

%% ! No IRS: R-E region by GP and SDR
% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');

% * GP and SDR
% ni_gp;
ni_sdr;
save('data/re_ni_gp_sdr.mat');

%% * R-E plots
figure('name', 'No IRS: R-E region by GP and SDR')
plot(niGpSample(1, :) / nSubbands, 1e6 * niGpSample(2, :), 'k--');
hold on;
plot(niSdrSample(1, :) / nSubbands, 1e6 * niSdrSample(2, :), 'r-');
hold off;
grid minor;
legend('GP', 'SDR');
xlabel('Rate per subband [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_ni_gp_sdr.fig');
