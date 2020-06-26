clear; clc; setup; config; load('data/tap.mat');

[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

ni_sdr;
ff_sdr;
fs_sdr;

%% * R-E plot
figure('name', 'R-E region by No IRS, FF-IRS and FS-IRS')
plot(niSdrSample(1, :), 1e6 * niSdrSample(2, :), 'k--');
hold on;
plot(ffSdrSample(1, :), 1e6 * ffSdrSample(2, :), 'r-');
hold on;
plot(fsSdrSample(1, :), 1e6 * fsSdrSample(2, :), 'b-.');
hold off;
grid minor;
legend('No IRS', 'FF-IRS', 'FS-IRS');
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
save('data/re_irs.mat');
savefig('plot/re_irs.fig');
