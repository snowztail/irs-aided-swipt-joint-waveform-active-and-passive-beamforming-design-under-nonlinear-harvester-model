clear; clc; setup; config; load('data/tap.mat');

% * Generate channels
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

% * Initialize algorithm
[capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);

rateConstraint = linspace(capacity, 0, nSamples);

ffSdrSolution = cell(nSamples, 1);
ffSdrSolution{1}.powerRatio = powerRatio;
ffSdrSolution{1}.infoWaveform = infoWaveform;
ffSdrSolution{1}.powerWaveform = powerWaveform;

ffSdrSample = zeros(2, nSamples);
ffSdrSample(:, 1) = [capacity; 0];

% * SDR
for iSample = 2 : nSamples
    % * AO
    isOuterConverged = false;
    current_ = 0;
    while ~isOuterConverged
        [irs] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint(iSample), nCandidates, tolerance);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        isInnerConverged = false;
        current__ = 0;
        while ~isInnerConverged
            [infoRatio, powerRatio] = split_ratio(compositeChannel, infoWaveform, noisePower, rateConstraint(iSample));
            [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), nCandidates, tolerance);
            [rate, current] = re_sample(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower);
            isInnerConverged = abs(current - current__) / current <= tolerance || current <= 1e-10;
            current__ = current;
        end
        isOuterConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end
    ffSdrSolution{iSample}.powerRatio = powerRatio;
    ffSdrSolution{iSample}.infoWaveform = infoWaveform;
    ffSdrSolution{iSample}.powerWaveform = powerWaveform;
    ffSdrSample(:, iSample) = [rate; current];
end
