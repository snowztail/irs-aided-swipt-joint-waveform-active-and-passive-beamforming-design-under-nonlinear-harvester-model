clear; clc; setup; config; load('data/tap.mat');

[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

%% ! IRS: R-E region by SDR
% * Initialize algorithm
[capacity, irs, infoWaveformOpt] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[current, ~, ~, powerWaveformOpt, infoRatio, powerRatio] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs, txPower, nCandidates, noisePower);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);

% * SDR
ffReSdrSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;

    % * Initialize waveform and splitting ratio for each sample
    infoWaveform = infoWaveformOpt;
    powerWaveform = powerWaveformOpt;
    while ~isConverged
        [infoRatio, powerRatio] = split_ratio(compositeChannel, noisePower, rateConstraint(iSample), infoWaveform);
        [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        [rate, current] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end

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
    ffReSdrSample(:, iSample) = [rate; current; powerRatio];
end

%% ! IRS: R-E region by GP
% * Initialize algorithm
[capacity, irs, infoWaveformOpt] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower);
[current, ~, ~, powerWaveformOpt, infoRatio, powerRatio] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs, txPower, nCandidates, noisePower);
[compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
rateConstraint = linspace(0, (1 - tolerance) * capacity, nSamples);

% * GP
ffReGpSample = zeros(3, nSamples);
for iSample = 1 : nSamples
    isConverged = false;
    current_ = 0;

    % * Initialize waveform and splitting ratio for each sample
    infoRatio = 1 - eps;
    powerRatio = eps;
    infoWaveform = infoWaveformOpt;
    powerWaveform = powerWaveformOpt;
    [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);

    % * AO
    while ~isConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint(iSample), tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint(iSample), tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end
    ffReGpSample(:, iSample) = [rate; current; powerRatio];
end

%% * R-E plot
figure('name', 'R-E region by GP and SDR')
plot(ffReGpSample(1, :), 1e6 * ffReGpSample(2, :), 'k--');
hold on;
plot(ffReSdrSample(1, :), 1e6 * ffReSdrSample(2, :), 'r-');
hold off;
grid minor;
legend('GP', 'SDR');
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
save('data/re_gp_sdr.mat');
savefig('plot/re_gp_sdr.fig');
