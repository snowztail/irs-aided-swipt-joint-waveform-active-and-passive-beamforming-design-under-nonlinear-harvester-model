clear; clc; setup; config_subband; load('data/tap.mat');

%% ! R-E region vs number of subbands
niSample = cell(2, length(Variable.nSubbands));
niSolution = cell(2, length(Variable.nSubbands));
ffSample = cell(2, length(Variable.nSubbands));
ffSolution = cell(2, length(Variable.nSubbands));

for iSubband = 1 : length(Variable.nSubbands)
    % * Update channels
    nSubbands = Variable.nSubbands(iSubband);
    [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

    % * GP and SDR
    ni_gp;
    niSample{1, iSubband} = niGpSample;
    niSolution{1, iSubband} = niGpSolution;
    ni_sdr;
    niSample{2, iSubband} = niSdrSample;
    niSolution{2, iSubband} = niSdrSolution;

    ff_gp;
    ffSample{1, iSubband} = ffGpSample;
    ffSolution{1, iSubband} = ffGpSolution;
    ff_sdr;
    ffSample{2, iSubband} = ffSdrSample;
    ffSolution{2, iSubband} = ffSdrSolution;
end
save('data/re_subband.mat');

%% * R-E plots
% * No IRS/
figure('name', 'No IRS: R-E region vs number of subbands');
legendString = cell(2 * length(Variable.nSubbands), 1);
% * GP
ax = gca;
ax.ColorOrderIndex = 1;
for iSubband = 1 : length(Variable.nSubbands)
    plot(niSample{1, iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * niSample{1, iSubband}(2, :));
    legendString{iSubband} = sprintf('GP: N = %d', Variable.nSubbands(iSubband));
    hold on;
end
% * SDR
ax = gca;
ax.ColorOrderIndex = 1;
for iSubband = 1 : length(Variable.nSubbands)
    plot(niSample{2, iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * niSample{2, iSubband}(2, :), '--');
    legendString{length(Variable.nSubbands) + iSubband} = sprintf('SDR: N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_subband_ni.fig');

% * IRS
figure('name', 'IRS: R-E region vs number of subbands');
legendString = cell(2 * length(Variable.nSubbands), 1);
% * GP
ax = gca;
ax.ColorOrderIndex = 1;
for iSubband = 1 : length(Variable.nSubbands)
    plot(ffSample{1, iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * ffSample{1, iSubband}(2, :));
    legendString{iSubband} = sprintf('GP: N = %d', Variable.nSubbands(iSubband));
    hold on;
end
% * SDR
ax = gca;
ax.ColorOrderIndex = 1;
for iSubband = 1 : length(Variable.nSubbands)
    plot(ffSample{2, iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * ffSample{2, iSubband}(2, :), '--');
    legendString{length(Variable.nSubbands) + iSubband} = sprintf('SDR: N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_subband_ff.fig');
