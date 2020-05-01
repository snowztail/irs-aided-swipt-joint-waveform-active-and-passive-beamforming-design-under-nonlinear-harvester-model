%% * Initialize script
clear; close all; clc; config_distances;

directRate = nan(1, length(rateConstraint));
compositeRate = nan(length(Variable.distance), length(rateConstraint));
directCurrent = nan(1, length(rateConstraint));
compositeCurrent = nan(length(Variable.distance), length(rateConstraint));

%% * Generate fading channels
[directFading] = channel_tgn_e(nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
[incidentFading] = channel_tgn_e(nReflectors, nSubbands, nUsers, carrierFrequency, fadingType);
reflectiveFading = zeros(nUsers, nSubbands * nReflectors);
for iReflector = 1 : nReflectors
    [reflectiveFading(iReflector : nReflectors : end)] = channel_tgn_e(nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
end
[directPathloss] = large_scale_fading(directDistance);
directChannel = directFading ./ sqrt(directPathloss);
directChannel = directChannel.';

%% * Waveform design without IRS
for iConstraint = 1 : length(rateConstraint)
    [~, ~, directRate(iConstraint), directCurrent(iConstraint)] = waveform_gp(k2, k4, resistance, txPower, noisePower, rateConstraint(iConstraint), tolerance, directChannel);
    if isnan(directCurrent(iConstraint))
        break;
    end
end

%% * Waveform design with IRS
for iDistance = 1 : length(Variable.distance)
    % distances
    incidentDistance = Variable.distance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    % pathlosses
    [incidentPathloss] = large_scale_fading(incidentDistance);
    [reflectivePathloss] = large_scale_fading(reflectiveDistance);
    % channels
    incidentChannel = incidentFading ./ sqrt(incidentPathloss);
    reflectiveChannel = reflectiveFading ./ sqrt(reflectivePathloss);
    [irs] = irs_selective(directChannel, incidentChannel, reflectiveChannel);
    [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    compositeChannel = compositeChannel.';
    % rate and current
    for iConstraint = 1 : length(rateConstraint)
        [~, ~, compositeRate(iDistance, iConstraint), compositeCurrent(iDistance, iConstraint)] = waveform_gp(k2, k4, resistance, txPower, noisePower, rateConstraint(iConstraint), tolerance, compositeChannel);
        if isnan(compositeCurrent(iDistance, iConstraint))
            break;
        end
    end
end

%% * Plot R-E region
legendString = cell(length(Variable.distance) + 1, 1);
figure('Name', 'R-E region vs AP-IRS distance');
plot(directRate, directCurrent * 1e6);
legendString{1} = sprintf('No IRS');
hold on;
for iDistance = 1 : length(Variable.distance)
    plot(compositeRate(iDistance, :), compositeCurrent(iDistance, :) * 1e6);
    legendString{iDistance + 1} = sprintf('IRS: d = %.2f', Variable.distance(iDistance));
    hold on;
end
hold off;
grid on; grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('I_{DC} [\muA]');
