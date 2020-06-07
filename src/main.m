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
while ~isConverged
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    [compositeCapacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(compositeChannel, subbandPower);
    [irs, maxRate] = irs_flat_wit(noisePower, concatMatrix, infoWaveform, nCandidates);
    isConverged = abs(maxRate - maxRate_) / maxRate <= tolerance;
    maxRate_ = maxRate;
end
nSamples = floor(maxRate);
rateConstraint = nSamples : -1 : 1;

% * Achievable R-E region by FF-IRS
rePair = zeros(2, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    rate_ = 0;
    while ~isConverged
        [irs] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint(iSample), tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);
        [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
        isConverged = (current - current_) / current <= tolerance || abs(rate - rate_) / rate <= tolerance;
        current_ = current;
        rate_ = rate;
    end
    rePair(:, iSample) = [current; rate];
end
