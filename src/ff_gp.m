%% ! FF-IRS: R-E region by GP
% * Initialize algorithm
[capacity, irs, infoWaveformOpt] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[current, ~, ~, powerWaveformOpt] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs, txPower, nCandidates, noisePower);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);

% * GP
ffGpSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;

    % * Initialize waveform and splitting ratio for each sample
    infoRatio = 1 - eps;
    powerRatio = eps;
    infoWaveform = infoWaveformOpt;
    powerWaveform = powerWaveformOpt;
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);

    % * AO
    while ~isConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end
    ffGpSample(:, iSample) = [rate; current; powerRatio];
end
