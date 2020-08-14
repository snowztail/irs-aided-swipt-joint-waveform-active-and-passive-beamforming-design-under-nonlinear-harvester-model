clear; clc; setup; config_los;

%% ! R-E region for IRS-aided NLOS and LOS channels
reSample = cell(nChannels, nCases);
reSolution = cell(nChannels, nCases);

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directNlosTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentNlosTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveNlosTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    [incidentLosTapGain, ~] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveLosTapGain, ~] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directNlosChannel] = frequency_response(directNlosTapGain, directTapDelay, directDistance, subbandFrequency, fadingMode);
    [incidentNlosChannel] = frequency_response(incidentNlosTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
    [reflectiveNlosChannel] = frequency_response(reflectiveNlosTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);

    [incidentLosChannel] = frequency_response(incidentLosTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
    [reflectiveLosChannel] = frequency_response(reflectiveLosTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);

    % * Optimization based on NLOS channels
    [reSample{iChannel, 1}, reSolution{iChannel, 1}] = re_sample(beta2, beta4, directNlosChannel, incidentNlosChannel, reflectiveNlosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

    % * Optimization based on NLOS direct channel and LOS incident and reflective channels
    [reSample{iChannel, 2}, reSolution{iChannel, 2}] = re_sample(beta2, beta4, directNlosChannel, incidentLosChannel, reflectiveLosChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end

% * Average over channel realizations
reSampleAvg = cell(1, nCases);
for iCase = 1 : nCases
    reSampleAvg{iCase} = mean(cat(3, reSample{:, iCase}), 3);
end
save('data/re_los.mat');

%% * R-E plots
figure('name', 'R-E region for IRS-aided NLOS and LOS channels');
for iCase = 1 : nCases
    plot(reSampleAvg{iCase}(1, :) / nSubbands, 1e6 * reSampleAvg{iCase}(2, :));
    hold on;
end
hold off;
grid minor;
legend('NLOS', 'LOS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_los.fig');
