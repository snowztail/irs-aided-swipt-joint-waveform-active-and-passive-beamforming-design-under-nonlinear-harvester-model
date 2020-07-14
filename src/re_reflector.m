clear; clc; setup; config_reflector; load('data/tap.mat');

%% ! R-E region vs number of IRS elements
% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), directDistance, directTapGain, directTapDelay, 'direct');
[incidentChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
[reflectiveChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');
gpSample = cell(length(Variable.nReflectors), 1);
sdrSample = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    % * Update channels
    nReflectors = Variable.nReflectors(iReflector);
    incidentChannel = incidentChannel_(:, :, 1 : nReflectors);
    reflectiveChannel = reflectiveChannel_(:, 1 : nReflectors, :);

    % * GP and SDR
    ff_gp;
    gpSample{iReflector} = ffGpSample;
    ff_sdr;
    sdrSample{iReflector} = ffSdrSample;
end
save('data/re_reflector.mat')

%% * R-E plots
figure('name', 'R-E region vs number of reflectors');
legendString = cell(2 * length(Variable.nReflectors), 1);

% * GP
ax = gca;
ax.ColorOrderIndex = 1;
for iReflector = 1 : length(Variable.nReflectors)
    plot(gpSample{iReflector}(1, :), 1e6 * gpSample{iReflector}(2, :));
    legendString{iReflector} = sprintf('GP: L = %d', Variable.nReflectors(iReflector));
    hold on;
end

% * SDR
ax = gca;
ax.ColorOrderIndex = 1;
for iReflector = 1 : length(Variable.nReflectors)
    plot(sdrSample{iReflector}(1, :), 1e6 * sdrSample{iReflector}(2, :), '--');
    legendString{length(Variable.nReflectors) + iReflector} = sprintf('SDR: L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_reflector.fig');
