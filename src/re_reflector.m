clear; clc; setup; config_reflector; load('data/tap.mat');

%% ! R-E region vs number of IRS elements
niSample = cell(2, 1);
niSolution = cell(2, 1);
ffSample = cell(2, length(Variable.nReflectors));
ffSolution = cell(2, length(Variable.nReflectors));

% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), directDistance, directTapGain, directTapDelay, 'direct');
[incidentChannelData] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
[reflectiveChannelData] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');

% * No IRS
ni_gp;
niSample{1} = niGpSample;
niSolution{1} = niGpSolution;
ni_sdr;
niSample{2} = niSdrSample;
niSolution{2} = niSdrSolution;

for iReflector = 1 : length(Variable.nReflectors)
    % * Update channels
    nReflectors = Variable.nReflectors(iReflector);
    incidentChannel = incidentChannelData(:, :, 1 : nReflectors);
    reflectiveChannel = reflectiveChannelData(:, 1 : nReflectors, :);

    % * GP and SDR
    ff_gp;
    ffSample{1, iReflector} = ffGpSample;
    ffSolution{1, iReflector} = ffGpSolution;
    ff_sdr;
    ffSample{2, iReflector} = ffSdrSample;
    ffSolution{2, iReflector} = ffSdrSolution;
end
save('data/re_reflector.mat')

%% * R-E plots
figure('name', 'R-E region vs number of reflectors');
legendString = cell(2 * length(Variable.nReflectors) + 2, 1);

% * GP
ax = gca;
ax.ColorOrderIndex = 1;
% * No IRS
plot(niSample{1}(1, :), 1e6 * niSample{1}(2, :));
legendString{1} = sprintf('GP: L = 0');
hold on;
for iReflector = 1 : length(Variable.nReflectors)
    plot(ffSample{1, iReflector}(1, :), 1e6 * ffSample{1, iReflector}(2, :));
    legendString{iReflector + 1} = sprintf('GP: L = %d', Variable.nReflectors(iReflector));
    hold on;
end

% * SDR
ax = gca;
ax.ColorOrderIndex = 1;
plot(niSample{2}(1, :), 1e6 * niSample{2}(2, :));
legendString{length(Variable.nReflectors) + 2} = sprintf('SDR: L = 0');
hold on;
for iReflector = 1 : length(Variable.nReflectors)
    plot(ffSample{2, iReflector}(1, :), 1e6 * ffSample{2, iReflector}(2, :), '--');
    legendString{length(Variable.nReflectors) + iReflector + 2} = sprintf('SDR: L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_reflector.fig');
