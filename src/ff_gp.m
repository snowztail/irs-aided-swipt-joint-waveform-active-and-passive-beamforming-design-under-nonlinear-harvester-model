%% ! FF-IRS: R-E region by GP
% * Initialize algorithm
[capacity, irs, infoWaveform, powerWaveform] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
infoRatio = 1 - eps;
powerRatio = 1 - infoRatio;

rateConstraint = linspace(capacity, 0, nSamples);

ffGpSolution = cell(nSamples, 1);
ffGpSolution{1}.powerRatio = powerRatio;
ffGpSolution{1}.infoWaveform = infoWaveform;
ffGpSolution{1}.powerWaveform = powerWaveform;

ffGpSample = zeros(2, nSamples);
ffGpSample(:, 1) = [capacity; 0];

% * GP
for iSample = 2 : nSamples
    isConverged = false;
    current_ = 0;
    % * AO
    while ~isConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end
    ffGpSolution{iSample}.powerRatio = powerRatio;
    ffGpSolution{iSample}.infoWaveform = infoWaveform;
    ffGpSolution{iSample}.powerWaveform = powerWaveform;
    ffGpSample(:, iSample) = [rate; current];
end
