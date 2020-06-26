clear; clc; setup; config; load('data/tap.mat');

[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

ff_gp;
ff_sdr;

%% * R-E plot
figure('name', 'FF-IRS: R-E region by GP and SDR')
plot(ffGpSample(1, :), 1e6 * ffGpSample(2, :), 'k--');
hold on;
plot(ffSdrSample(1, :), 1e6 * ffSdrSample(2, :), 'r-');
hold off;
grid minor;
legend('GP', 'SDR');
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
save('data/re_ff_gp_sdr.mat');
savefig('plot/re_ff_gp_sdr.fig');
