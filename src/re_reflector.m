clear; clc; setup; config_reflector;

%% ! R-E region vs number of IRS elements
reSample = cell(length(Variable.nReflectors), 1);
reSolution = cell(length(Variable.nReflectors), 1);

for iReflector = 1 : length(Variable.nReflectors)
    % * Generate channels
    nReflectors = Variable.nReflectors(iReflector);
    [directChannel] = frequency_response(Variable.directTapGain{iReflector}, directTapDelay, directDistance, subbandFrequency, fadingMode);
    [incidentChannel] = frequency_response(Variable.incidentTapGain{iReflector}, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
    [reflectiveChannel] = frequency_response(Variable.reflectiveTapGain{iReflector}, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);

    % * Alternating optimization
    [reSample{iReflector}, reSolution{iReflector}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
end
save('data/re_reflector.mat');

%% * R-E plots
figure('name', 'R-E region vs number of reflectors');
legendString = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    plot(reSample{iReflector}(1, :) / nSubbands, 1e6 * reSample{iReflector}(2, :));
    legendString{iReflector} = sprintf('L = %d', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('plots/re_reflector.fig');
