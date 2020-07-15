%% ! No IRS: R-E region by SDR
% * Initialize algorithm
[capacity, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_fs(directChannel, txPower, noisePower);

rateConstraint = linspace(capacity, 0, nSamples);

niSdrSolution = cell(nSamples, 1);
niSdrSolution{1}.powerRatio = powerRatio;
niSdrSolution{1}.infoWaveform = infoWaveform;
niSdrSolution{1}.powerWaveform = powerWaveform;

niSdrSample = zeros(2, nSamples);
niSdrSample(:, 1) = [capacity; 0];

% * SDR
for iSample = 2 : nSamples
    % * AO
    isConverged = false;
    current_ = 0;
    while ~isConverged
        [infoRatio, powerRatio] = split_ratio(directChannel, noisePower, rateConstraint(iSample), infoWaveform);
        [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, directChannel, infoWaveform, powerWaveform);
        [rate, current] = re_sample(beta2, beta4, directChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end
    niSdrSolution{iSample}.powerRatio = powerRatio;
    niSdrSolution{iSample}.infoWaveform = infoWaveform;
    niSdrSolution{iSample}.powerWaveform = powerWaveform;
    niSdrSample(:, iSample) = [rate; current];
end
