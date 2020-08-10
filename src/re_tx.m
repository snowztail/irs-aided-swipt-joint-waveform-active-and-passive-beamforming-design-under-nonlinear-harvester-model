clear; clc; setup; config_tx;

%% ! R-E region vs number of transmit antennas
reSample = cell(length(Variable.nTxs), 1);
reSolution = cell(length(Variable.nTxs), 1);

for iTx = 1 : length(Variable.nTxs)
    % * Generate channels
    [directChannel] = frequency_response(Variable.directTapGain{iTx}, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
    [incidentChannel] = frequency_response(Variable.incidentTapGain{iTx}, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
    [reflectiveChannel] = frequency_response(Variable.reflectiveTapGain{iTx}, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

    % * Alternating optimization
    [reSample{iTx}, reSolution{iTx}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end
save('data/re_tx.mat');

%% * R-E plots
figure('name', 'R-E region vs number of transmit antennas');
legendString = cell(length(Variable.nTxs), 1);
for iTx = 1 : length(Variable.nTxs)
    plot(reSample{iTx}(1, :) / nSubbands, 1e6 * reSample{iTx}(2, :));
    legendString{iTx} = sprintf('M = %d', Variable.nTxs(iTx));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_tx.fig');
