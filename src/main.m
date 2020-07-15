clear; clc; setup; config; load('data/tap.mat');

% * Generate channels
[directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, nReflectors, subbandFrequency, fadingMode, 'direct');
[incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, nReflectors, subbandFrequency, fadingMode, 'incident');
[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, nReflectors, subbandFrequency, fadingMode, 'reflective');

[capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
irs = exp(1i * 2 * pi * rand(nReflectors, 1));

[compositeChannel, concatChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

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
    ffSdrSolution{iSample}.powerRatio = powerRatio;
    ffSdrSolution{iSample}.infoWaveform = infoWaveform;
    ffSdrSolution{iSample}.powerWaveform = powerWaveform;
    ffSdrSample(:, iSample) = [rate; current];
end
