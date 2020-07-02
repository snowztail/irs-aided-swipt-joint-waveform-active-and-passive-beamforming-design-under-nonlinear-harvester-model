%% ! FF-IRS: R-E region by SDR
% * Initialize algorithm
[capacity, irs, infoWaveformOpt] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[current, ~, ~, powerWaveformOpt] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs, txPower, nCandidates, noisePower);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);

% * SDR
ffSdrSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    % * Initialize waveform and splitting ratio for each sample
    infoWaveform = infoWaveformOpt;
    powerWaveform = powerWaveformOpt;

    isConverged = false;
    current_ = 0;
    while ~isConverged
        [infoRatio, powerRatio] = split_ratio(compositeChannel, noisePower, rateConstraint(iSample), infoWaveform);
        [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        [rate, current] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end

    % * AO
    isOuterConverged = false;
    current_ = 0;
    while ~isOuterConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        isInnerConverged = false;
        current__ = 0;
        while ~isInnerConverged
            [infoRatio, powerRatio] = split_ratio(compositeChannel, noisePower, rateConstraint(iSample), infoWaveform);
            [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
            [rate, current] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
            isInnerConverged = abs(current - current__) / current <= tolerance || current <= 1e-10;
            current__ = current;
        end
        isOuterConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end
    ffSdrSample(:, iSample) = [rate; current; powerRatio];
end