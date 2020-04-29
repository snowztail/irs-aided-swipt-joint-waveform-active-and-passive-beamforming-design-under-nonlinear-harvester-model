%% * Initialize script
clear; close all; clc; config_reflectors;

directRate = nan(1, length(rateConstraint));
compositeRate = nan(length(Variable.nReflectors), length(rateConstraint));
directCurrent = nan(1, length(rateConstraint));
compositeCurrent = nan(length(Variable.nReflectors), length(rateConstraint));

%% * Generate channels
% K * N
[directFading] = channel_tgn_e(nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
% L * N
[incidentFading] = channel_tgn_e(Variable.nReflectors(end), nSubbands, nUsers, carrierFrequency, fadingType);
% K * NL
reflectiveFading = zeros(nUsers, nSubbands * Variable.nReflectors(end));
for iReflector = 1 : Variable.nReflectors(end)
    [reflectiveFading(iReflector : Variable.nReflectors(end) : end)] = channel_tgn_e(nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
end
directChannel = directFading ./ sqrt(directPathloss);

%% * Waveform design without IRS
for iConstraint = 1 : length(rateConstraint)
    [~, ~, directRate(iConstraint), directCurrent(iConstraint)] = waveform_gp(k2, k4, resistance, txPower, noisePower, rateConstraint(iConstraint), tolerance, directChannel);
end

%% * Waveform design with IRS
for iReflector = 1 : length(Variable.nReflectors)
    incidentChannel = incidentFading(1 : Variable.nReflectors(iReflector), :) ./ sqrt(incidentPathloss);
    reflectiveChannel = zeros(Variable.nReflectors(iReflector), nSubbands);
    for jReflector = 1 : Variable.nReflectors(iReflector)
        for iSubband = 1 : nSubbands
            reflectiveChannel(jReflector + (iSubband - 1) * Variable.nReflectors(iReflector)) = reflectiveFading(jReflector + (iSubband - 1) * Variable.nReflectors(end)) ./ sqrt(incidentPathloss);
        end
    end
    [irs] = irs_selective(directChannel, incidentChannel, reflectiveChannel);
    [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    compositeChannel = compositeChannel.';
    % rate and current
    for iConstraint = 1 : length(rateConstraint)
        [~, ~, compositeRate(iReflector, iConstraint), compositeCurrent(iReflector, iConstraint)] = waveform_gp(k2, k4, resistance, txPower, noisePower, rateConstraint(iConstraint), tolerance, compositeChannel);
        if isnan(compositeCurrent(iReflector, iConstraint))
            break;
        end
    end
end

%% * Plot R-E region
legendString = cell{length(Variable.nReflectors) + 1, 1};
figure('Name', 'R-E region vs number of reflectors');
plot(directRate, directCurrent * 1e6);
legendString{1} = sprintf('Plain');
hold on;
for iReflector = 1 : length(Variable.nReflectors)
    plot(compositeRate(iReflector, :), compositeCurrent(iReflector, :) * 1e6);
    legendString{iReflector + 1} = sprintf('IRS: L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid on; grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('I_{DC} [\muA]');
ylim([0, 6]);
