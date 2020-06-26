clear; clc; setup; config_subband; load('data/tap.mat');

%% ! R-E region vs number of subbands
niSample = cell(length(Variable.nSubbands), 1);
ffSample = cell(length(Variable.nSubbands), 1);
fsSample = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    % * Update channels
    nSubbands = Variable.nSubbands(iSubband);
    [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

    % * SDR
    ni_sdr;
    niSample{iSubband} = niSdrSample;
    ff_sdr;
    ffSample{iSubband} = ffSdrSample;
    fs_sdr;
    fsSample{iSubband} = fsSdrSample;
end
save('data/re_subband.mat');

%% * R-E plots
figure('name', 'No IRS: R-E region vs number of subbands');
legendString = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    plot(niSample{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * niSample{iSubband}(2, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_subband_ni.fig');

figure('name', 'FF-IRS: R-E region vs number of subbands');
legendString = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    plot(ffSample{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * ffSample{iSubband}(2, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_subband_ff.fig');

figure('name', 'FS-IRS: R-E region vs number of subbands');
legendString = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    plot(fsSample{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * fsSample{iSubband}(2, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_subband_fs.fig');
