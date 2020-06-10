clear; clc; setup; config;

% * Direct link
[directTapGain, directTapDelay] = taps_tgn(nTxs, nRxs);
[directFading] = fading_tgn(directTapGain, directTapDelay, nSubbands, subbandFrequency, fadingMode);
[directPathloss] = path_loss(directDistance, "direct");
directChannel = directFading / sqrt(directPathloss);

% * Incident link
[incidentTapGain, incidentTapDelay] = taps_tgn(nTxs, nReflectors);
[incidentFading] = fading_tgn(incidentTapGain, incidentTapDelay, nSubbands, subbandFrequency, fadingMode);
[incidentPathloss] = path_loss(incidentDistance, "incident");
incidentChannel = incidentFading / sqrt(incidentPathloss);

% * Reflective link
[reflectiveTapGain, reflectiveTapDelay] = taps_tgn(nReflectors, nRxs);
[reflectiveFading] = fading_tgn(reflectiveTapGain, reflectiveTapDelay, nSubbands, subbandFrequency, fadingMode);
[reflectivePathloss] = path_loss(reflectiveDistance, "reflective");
reflectiveChannel = reflectiveFading / sqrt(reflectivePathloss);

% * Initialize algorithm by WIT
isConverged = false;
maxRate_ = 0;
irs = irsGain * ones(nReflectors, 1);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
[compositeCapacity_, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
[infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(compositeChannel, subbandPower);
while ~isConverged
    [irs, maxRate] = irs_flat_wit(noisePower, concatMatrix, infoWaveform, nCandidates);
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    [compositeCapacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(compositeChannel, subbandPower);
    isConverged = abs(compositeCapacity - compositeCapacity_) / compositeCapacity <= tolerance;
    compositeCapacity_ = compositeCapacity;
end
rateConstraint = resolution * (floor(maxRate / resolution) : -1 : 0);
nSamples = length(rateConstraint);

% * Achievable R-E region by FF-IRS
reSample = zeros(2, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    rate_ = 0;
    [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
    while ~isConverged
        [irs, currentInit, rateInit] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint(iSample), tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);
        [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
        isConverged = (current - current_) / current <= tolerance || abs(rate - rate_) / rate <= tolerance;
        current_ = current;
        rate_ = rate;
    end
    reSample(:, iSample) = [current; rate];
end
flag = 1;
