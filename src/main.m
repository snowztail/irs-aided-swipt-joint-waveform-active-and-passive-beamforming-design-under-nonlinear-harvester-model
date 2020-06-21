clear; clc; setup; config; load('data/tap.mat');

[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

%% ! No-IRS: R-E region
% * Initialize algorithm by WIT
[capacity, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_no_irs(directChannel, txPower, noisePower);
rateConstraint = capacity : -capacity / (nSamples - 1) : 0;

% * Achievable R-E region without IRS
directReSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    while ~isConverged
        [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, directChannel, infoWaveform, powerWaveform);
        [infoRatio, powerRatio] = split_ratio(infoWaveform, noisePower, rateConstraint(iSample), directChannel);
        [rate, current] = re_sample(beta2, beta4, directChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
        isConverged = abs(current - current_) / current <= tolerance || current == 0;
        current_ = current;
    end
    directReSample(:, iSample) = [rate; current; powerRatio];
end

%% ! IRS: R-E region
% * Initialize algorithm by WIT
[capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
rateConstraint = capacity : -capacity / (nSamples - 1) : 0;
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
% * Achievable R-E region by FF-IRS
ffReSample = zeros(3, nSamples);
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
    ffReSample(:, iSample) = [rate; current; powerRatio];
end
