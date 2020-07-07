clear; clc; setup; config_irs; load('data/tap.mat');

%% ! R-E region by no IRS, FF-IRS and FS-IRS
% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

% * SDR
ni_sdr;
ff_sdr;
fs_sdr;
save('data/re_irs.mat');

%% * R-E plots
figure('name', 'R-E region by no IRS, FF-IRS and FS-IRS');
plot(niSdrSample(1, :) / nSubbands, 1e6 * niSdrSample(2, :), 'k--');
hold on;
plot(ffSdrSample(1, :) / nSubbands, 1e6 * ffSdrSample(2, :), 'r-');
hold on;
plot(fsSdrSample(1, :) / nSubbands, 1e6 * fsSdrSample(2, :), 'b-.');
hold off;
grid minor;
legend('No IRS', 'FF-IRS', 'FS-IRS');
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_irs.fig');
