clear; clc; setup; config_subband;

%% ! R-E region vs number of subbands
reSample = cell(nChannels, length(Variable.nSubbands));
reSolution = cell(nChannels, length(Variable.nSubbands));

for iChannel = 1 : nChannels
    for iSubband = 1 : length(Variable.nSubbands)
        % * Get number of subbands and subband frequency
        nSubbands = Variable.nSubbands(iSubband);
        [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);

        % * Generate tap gains and delays
        [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
        [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
        [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

        % * Construct channels
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

        % * Alternating optimization
        [reSample{iChannel, iSubband}, reSolution{iChannel, iSubband}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reSampleAvg = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
    reSampleAvg{iSubband} = mean(cat(3, reSample{:, iSubband}), 3);
end
save('data/re_subband.mat');

%% * R-E plots
figure('name', 'R-E region vs number of subbands');
legendString = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    plot(reSampleAvg{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reSampleAvg{iSubband}(2, :));
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
