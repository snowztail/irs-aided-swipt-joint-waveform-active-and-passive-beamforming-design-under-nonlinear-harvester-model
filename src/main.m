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
[infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(compositeChannel, subbandPower);
while ~isConverged
    [irs, maxRate] = irs_flat_wit(noisePower, concatMatrix, infoWaveform, nCandidates);
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    [compositeCapacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(compositeChannel, subbandPower);
    isConverged = abs(compositeCapacity - compositeCapacity_) / compositeCapacity <= tolerance;
    compositeCapacity_ = compositeCapacity;
end
rateConstraint = resolution * (floor(maxRate / resolution) : -1 : 0);
% rateConstraint = resolution * (0 : floor(maxRate / resolution));
nSamples = length(rateConstraint);

% % * Initialize algorithm by WPT
% [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wpt(compositeChannel, txPower);

% * Achievable R-E region by FF-IRS
reSample = zeros(2, nSamples);
infoRatio = 0.9; powerRatio = 0.1;
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;
    rate_ = 0;
%     powerRatio = iSample / nSamples; infoRatio = 1 - powerRatio;
%     [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wpt(compositeChannel, txPower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
    while ~isConverged
        [irs, currentInit, rateInit] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint(iSample), tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);
        [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
        % isConverged = abs(current - current_) / current <= tolerance && abs(rate - rate_) / rate <= tolerance;
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
        rate_ = rate;
    end
    reSample(:, iSample) = [current; rate];
%     [compositeCapacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
%     [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wpt(compositeChannel, txPower);
%     powerRatio = (iSample - 1) / nSamples; infoRatio = 1 - powerRatio;
end
flag = 1;
