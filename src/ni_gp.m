%% ! FF-IRS: R-E region by GP
% * Initialize algorithm
[capacity, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_fs(directChannel, txPower, noisePower);

rateConstraint = linspace(capacity, 0, nSamples);

niGpSolution = cell(nSamples, 1);
niGpSolution{1}.powerRatio = powerRatio;
niGpSolution{1}.infoWaveform = infoWaveform;
niGpSolution{1}.powerWaveform = powerWaveform;

niGpSample = zeros(2, nSamples);
niGpSample(:, 1) = [capacity; 0];

% * GP
for iSample = 2 : nSamples
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, directChannel, infoWaveform, powerWaveform);
    niGpSolution{iSample}.powerRatio = powerRatio;
    niGpSolution{iSample}.infoWaveform = infoWaveform;
    niGpSolution{iSample}.powerWaveform = powerWaveform;
    niGpSample(:, iSample) = [rate; current];
end
