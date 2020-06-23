clear; clc; setup; config_distance; load('data/tap.mat');

%% ! IRS: R-E region vs AP-IRS distance
ffReSample = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    % * Update channels
    incidentDistance = Variable.incidentDistance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

    % * Initialize algorithm by WIT
    [capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
    rateConstraint = capacity : -capacity / (nSamples - 1) : 0;
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Achievable R-E region by FF-IRS
    ffReSample{iDistance} = zeros(3, nSamples);
    for iSample = 1 : nSamples
        isConverged = false;
        current_ = 0;
        while ~isConverged
            [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
            [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
            [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
            [infoRatio, powerRatio] = split_ratio(infoWaveform, noisePower, rateConstraint(iSample), compositeChannel);
            [rate, current] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
            isConverged = abs(current - current_) / current <= tolerance || current == 0;
            current_ = current;
        end
        ffReSample{iDistance}(:, iSample) = [rate; current; powerRatio];
    end
end
