%% ! No IRS: R-E region by SDR
% * Initialize algorithm
[capacity, infoWaveform_] = wit_fs(directChannel, txPower, noisePower);
[current, ~, powerWaveform_] = wpt_fs(beta2, beta4, tolerance, directChannel, txPower, nCandidates, noisePower);
rateConstraint = linspace((1 - tolerance) * capacity, 0, nSamples);
% rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);
    infoWaveform = infoWaveform_;
    powerWaveform = zeros(size(powerWaveform_));
%     infoWaveform = infoWaveform_;
%     powerWaveform = powerWaveform_;

% * SDR
niSdrSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    % * Initialize waveform and splitting ratio for each sample
%     infoRatio = 1 - eps;
%     powerRatio = eps;
%     infoWaveform = infoWaveform_;
%     powerWaveform = powerWaveform_;
%     powerWaveform = zeros(size(powerWaveform_));
    
    [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, directChannel, infoWaveform, powerWaveform);
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
    niSdrSample(:, iSample) = [rate; current; powerRatio];
end
