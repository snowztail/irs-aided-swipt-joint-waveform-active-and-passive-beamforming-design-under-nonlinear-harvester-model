clear; clc; setup; config_subband;

%% ! R-E region vs number of subbands
reSample = cell(length(Variable.nSubbands), 1);
reSolution = cell(length(Variable.nSubbands), 1);

for iSubband = 1 : length(Variable.nSubbands)
    % * Generate channels
    nSubbands = Variable.nSubbands(iSubband);
    [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
    [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

    % * Alternating optimization
    [reSample{iSubband}, reSolution{iSubband}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end
save('data/re_subband.mat');

%% * R-E plots
figure('name', 'R-E region vs number of subbands');
legendString = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    plot(reSample{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reSample{iSubband}(2, :));
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_subband.fig');
