clear; clc; setup; config_reflector; load('data/tap.mat');

%% ! R-E region vs number of IRS elements
% * Generate channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), directDistance, directTapGain, directTapDelay, 'direct');
[incidentChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), incidentDistance, incidentTapGain, incidentTapDelay, 'incident');
[reflectiveChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, max(Variable.nReflectors), reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, 'reflective');
ffSample = cell(length(Variable.nReflectors), 1);
fsSample = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    % * Update channels
    nReflectors = Variable.nReflectors(iReflector);
    incidentChannel = incidentChannel_(:, :, 1 : nReflectors);
    reflectiveChannel = reflectiveChannel_(:, 1 : nReflectors, :);

    % * SDR
    ff_sdr;
    ffSample{iReflector} = ffSdrSample;
    fs_sdr;
    fsSample{iReflector} = fsSdrSample;
end
save('data/re_reflector.mat')

%% * R-E plots
figure('name', 'FF-IRS: R-E region vs number of reflectors');
legendString = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    plot(ffSample{iReflector}(1, :), 1e6 * ffSample{iReflector}(2, :));
    legendString{iReflector} = sprintf('L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_reflector_ff.fig');

figure('name', 'FS-IRS: R-E region vs number of reflectors');
legendString = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    plot(fsSample{iReflector}(1, :), 1e6 * fsSample{iReflector}(2, :));
    legendString{iReflector} = sprintf('L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
savefig('plots/re_reflector_fs.fig');
