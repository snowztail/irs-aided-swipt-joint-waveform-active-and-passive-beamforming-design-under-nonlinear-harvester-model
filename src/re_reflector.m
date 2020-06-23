clear; clc; setup; config_reflector; load('data/tap.mat');

% * Construct channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

%% ! IRS: R-E region vs number of IRS elements
ffReSample = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    % * Update channels
    nReflectors = Variable.nReflectors(iReflector);
    incidentChannel = incidentChannel_(:, :, 1 : nReflectors);
    reflectiveChannel = reflectiveChannel_(:, 1 : nReflectors, :);

    % * Initialize algorithm by WIT
    [capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
    rateConstraint = capacity : -capacity / (nSamples - 1) : 0;
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Achievable R-E region by FF-IRS
    ffReSample{iReflector} = zeros(3, nSamples);
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
        ffReSample{iReflector}(:, iSample) = [rate; current; powerRatio];
    end
end
