clear; clc; setup; config_distance; load('data/taps.mat');

%% ! R-E region vs AP-IRS distance
reSample = cell(length(Variable.incidentDistance), 1);
reSolution = cell(length(Variable.incidentDistance), 1);

% * Get tap data
directTapGain = directTapGain(:, 1 : nTxs);
incidentTapGain = incidentTapGain(:, 1 : nTxs, :);

% * Generate direct channel
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');

for iDistance = 1 : length(Variable.incidentDistance)
    % * Generate extra channels
    incidentDistance = Variable.incidentDistance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
    [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

    % * Alternating optimization
    [reSample{iDistance}, reSolution{iDistance}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end
save('data/re_distance.mat');

%% * R-E plots
figure('name', 'R-E region vs AP-IRS distance');
legendString = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    plot(reSample{iDistance}(1, :) / nSubbands, 1e6 * reSample{iDistance}(2, :));
    legendString{iDistance} = sprintf('d_I = %d', Variable.incidentDistance(iDistance));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_distance.fig');
