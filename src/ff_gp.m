%% ! FF-IRS: R-E region by GP
% * Initialize algorithm
[capacity, irs_, infoWaveform_] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[current, irs__, ~, powerWaveform_] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs_, txPower, nCandidates, noisePower);
[compositeChannel_, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs_);
rateConstraint = linspace((1 - tolerance) * capacity, 0, nSamples);
    infoRatio = 1 - eps;
    powerRatio = 1 - infoRatio;
    irs = irs_;
    compositeChannel = compositeChannel_;
    infoWaveform = infoWaveform_;
    powerWaveform = powerWaveform_;


% * GP
ffGpSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    % * Initialize waveform and splitting ratio for each sample
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);

%     infoRatio = 1 - eps;
%     powerRatio = 1 - infoRatio;

    isConverged = false;
    current_ = 0;
    % * AO
    while ~isConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
%         isConverged = true;
    end
    ffGpSample(:, iSample) = [rate; current; powerRatio];
end
