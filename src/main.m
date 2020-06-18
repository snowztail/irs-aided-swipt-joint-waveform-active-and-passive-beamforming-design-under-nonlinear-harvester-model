clear; clc; setup; config; load('data/tap.mat');

[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

%% ! No-IRS: R-E region
% % * Initialize algorithm by WIT
% [directCapacity, subbandPower] = channel_capacity(directChannel, txPower, noisePower);
% [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(directChannel, subbandPower);
% rateConstraint = directCapacity : -directCapacity / (nSamples - 1) : 0;
%
% % * Achievable R-E region without IRS
% directReSample = zeros(2, nSamples);
% for iSample = 1 : nSamples
%     [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, directChannel, nCandidates, scaleFactor);
%     directReSample(:, iSample) = [current; rate];
% end

%% ! IRS: R-E region
% * Initialize algorithm by WIT
isConverged = false;
maxRate_ = 0;
irs = irsGain * ones(nReflectors, 1);
while ~isConverged
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    [compositeCapacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(compositeChannel, subbandPower);
    [irs, maxRate] = irs_flat_wit(noisePower, concatMatrix, infoWaveform, nCandidates);
    isConverged = abs(maxRate - maxRate_) / maxRate <= tolerance;
    maxRate_ = maxRate;
end
rateConstraint = maxRate : -maxRate / (nSamples - 1) : 0;

% * Achievable R-E region by FF-IRS
ffReSample = zeros(2, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates, scaleFactor);
    while ~isConverged
        [irs, currentInit, rateInit] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint(iSample), tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);
        [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates, scaleFactor);
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end
    ffReSample(:, iSample) = [current; rate];
end
