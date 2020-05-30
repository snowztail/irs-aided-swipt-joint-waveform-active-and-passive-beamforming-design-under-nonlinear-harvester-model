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

% * Composite channel
irs = irsGain * ones(nReflectors, 1);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

% * Achievable rate with direct link only
[directRate] = channel_capacity(directChannel, txPower, noisePower);

% * Initialize algorithm
[infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(txPower, compositeChannel);

% * Achievable rate by FF-IRS
[irsMatrix, flatRate] = irs_sdr(beta2, beta4, noisePower, rateConstraint, tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates);

% * Achievable rate by FS-IRS
[irs_, compositeChannel_] = irs_selective(directChannel, incidentChannel, reflectiveChannel);
[selectiveRate] = channel_capacity(compositeChannel_, txPower, noisePower);
