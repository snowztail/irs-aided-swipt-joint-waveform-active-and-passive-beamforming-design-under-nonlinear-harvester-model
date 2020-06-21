clear; clc; setup; config_subband;

% * Tap response
[directTapGain, directTapDelay] = tap_tgn(nTxs, nRxs);
[incidentTapGain, incidentTapDelay] = tap_tgn(nTxs, nReflectors);
[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(nReflectors, nRxs);

% * R-E region vs number of subbands
directReSample = cell(length(Variable.nSubbands), 1);
ffReSample = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    nSubbands = Variable.nSubbands(iSubband);
    [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);

    % * Direct link
    [directFading] = fading_tgn(directTapGain, directTapDelay, nSubbands, subbandFrequency, fadingMode);
    [directPathloss] = path_loss(directDistance, "direct");
    directChannel = directFading / sqrt(directPathloss);

    % * Incident link
    [incidentFading] = fading_tgn(incidentTapGain, incidentTapDelay, nSubbands, subbandFrequency, fadingMode);
    [incidentPathloss] = path_loss(incidentDistance, "incident");
    incidentChannel = incidentFading / sqrt(incidentPathloss);

    % * Reflective link
    [reflectiveFading] = fading_tgn(reflectiveTapGain, reflectiveTapDelay, nSubbands, subbandFrequency, fadingMode);
    [reflectivePathloss] = path_loss(reflectiveDistance, "reflective");
    reflectiveChannel = reflectiveFading / sqrt(reflectivePathloss);

    %% ! No-IRS: R-E region vs number of subbands
    % * Initialize algorithm by WIT
    [directCapacity, subbandPower] = channel_capacity(directChannel, txPower, noisePower);
    [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(directChannel, subbandPower);
    rateConstraint = directCapacity : -directCapacity / (nSamples - 1) : 0;

    % * Achievable R-E region without IRS
    directReSample{iSubband} = zeros(2, nSamples);
    for iSample = 1 : nSamples
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, directChannel, nCandidates);
        directReSample{iSubband}(:, iSample) = [current; rate];
    end

    %% ! IRS: R-E region vs number of subbands
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
    ffReSample{iSubband} = zeros(2, nSamples);
    for iSample = 1 : nSamples
        isConverged = false;
        current_ = 0;
        [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
        while ~isConverged
            [irs, currentInit, rateInit] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint(iSample), tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);
            [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
            [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint(iSample), tolerance, compositeChannel, nCandidates);
            isConverged = abs(current - current_) / current <= tolerance || current == 0;
            current_ = current;
        end
        ffReSample{iSubband}(:, iSample) = [current; rate];
    end
end
