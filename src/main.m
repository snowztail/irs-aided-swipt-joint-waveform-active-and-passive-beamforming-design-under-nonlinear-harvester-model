clear; clc; setup; config; load('data/taps.mat');

% * Get tap data
directTapGain = directTapGain(:, 1 : nTxs);
incidentTapGain = incidentTapGain(:, 1 : nTxs, :);

% * Generate channels
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

% * Initialize algorithm and set rate constraints
[capacity, irs] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);

% * R-E sample
solution = cell(nSamples, 1);
sample = zeros(2, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(compositeChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_gp(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), tolerance);
    while ~isConverged
        [irs] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint(iSample), nCandidates, tolerance);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_gp(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), tolerance);
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
    end
    sample(:, iSample) = [rate; current];
    solution{iSample}.infoWaveform = infoWaveform;
    solution{iSample}.powerWaveform = powerWaveform;
    solution{iSample}.powerRatio = powerRatio;
end
