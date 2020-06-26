clear; clc; setup; config_distance; load('data/tap.mat');

%% ! R-E region vs AP-IRS distance
ffSample = cell(length(Variable.incidentDistance), 1);
fsSample = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    % * Update channels
    incidentDistance = Variable.incidentDistance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

    % * SDR
    ff_sdr;
    ffSample{iDistance} = ffSdrSample;
    fs_sdr;
    fsSample{iDistance} = fsSdrSample;
end
save('data/re_distance.mat');

%% * R-E plots
figure('name', 'FF-IRS: R-E region vs AP-IRS distance');
legendString = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    plot(ffSample{iDistance}(1, :), 1e6 * ffSample{iDistance}(2, :));
    legendString{iDistance} = sprintf('d_I = %d', Variable.incidentDistance(iDistance));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_distance_ff.fig');

figure('name', 'FS-IRS: R-E region vs AP-IRS distance');
legendString = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    plot(fsSample{iDistance}(1, :), 1e6 * fsSample{iDistance}(2, :));
    legendString{iDistance} = sprintf('d_I = %d', Variable.incidentDistance(iDistance));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_distance_fs.fig');
