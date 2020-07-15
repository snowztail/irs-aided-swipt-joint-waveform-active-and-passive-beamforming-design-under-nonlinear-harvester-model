clear; clc; setup; config_distance; load('data/tap.mat');

%% ! R-E region vs AP-IRS distance
niSample = cell(2, 1);
niSolution = cell(2, 1);
ffSample = cell(2, length(Variable.incidentDistance));
ffSolution = cell(2, length(Variable.incidentDistance));

% * Generate direct channel
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, 'direct');

% * No IRS
ni_gp;
niSample{1} = niGpSample;
niSolution{1} = niGpSolution;
ni_sdr;
niSample{2} = niSdrSample;
niSolution{2} = niSdrSolution;

for iDistance = 1 : length(Variable.incidentDistance)
    % * Update IRS-aided channels
    incidentDistance = Variable.incidentDistance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

    % * GP and SDR
    ff_gp;
    ffSample{1, iDistance} = ffGpSample;
    ffSolution{1, iDistance} = ffGpSolution;
    ff_sdr;
    ffSample{2, iDistance} = ffSdrSample;
    ffSolution{2, iDistance} = ffSdrSolution;
end
save('data/re_distance.mat');

%% * R-E plots
figure('name', 'R-E region vs AP-IRS distance');
legendString = cell(2 * length(Variable.incidentDistance) + 2, 1);

% * GP
ax = gca;
ax.ColorOrderIndex = 1;
% * No IRS
plot(niSample{1}(1, :), 1e6 * niSample{1}(2, :));
legendString{1} = sprintf('GP: No IRS');
hold on;
for iDistance = 1 : length(Variable.incidentDistance)
    plot(ffSample{1, iDistance}(1, :), 1e6 * ffSample{1, iDistance}(2, :));
    legendString{iDistance + 1} = sprintf('GP: d_I = %d', Variable.incidentDistance(iDistance));
    hold on;
end

% * SDR
ax = gca;
ax.ColorOrderIndex = 1;
plot(niSample{2}(1, :), 1e6 * niSample{2}(2, :));
legendString{length(Variable.incidentDistance) + 2} = sprintf('SDR: No IRS');
hold on;
for iDistance = 1 : length(Variable.incidentDistance)
    plot(ffSample{2, iDistance}(1, :), 1e6 * ffSample{2, iDistance}(2, :), '--');
    legendString{length(Variable.incidentDistance) + iDistance + 2} = sprintf('SDR: d_I = %d', Variable.incidentDistance(iDistance));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_distance.fig');
