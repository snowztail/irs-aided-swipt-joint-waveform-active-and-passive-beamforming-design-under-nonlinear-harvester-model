%% ! FF-IRS: R-E region by GP
% * Initialize algorithm
[capacity, infoWaveform_] = wit_fs(directChannel, txPower, noisePower);
[current, ~, powerWaveform_] = wpt_fs(beta2, beta4, tolerance, directChannel, txPower, nCandidates, noisePower);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);
    infoRatio = 1 - eps;
    powerRatio = eps;
    infoWaveform = infoWaveform_;
    powerWaveform = powerWaveform_;

% * GP
niGpSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    % * Initialize waveform and splitting ratio for each sample
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, directChannel, infoWaveform, powerWaveform);
    niGpSample(:, iSample) = [rate; current; powerRatio];
end
