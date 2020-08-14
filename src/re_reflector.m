clear; clc; setup; config_reflector;

%% ! R-E region vs number of IRS elements
reSample = cell(nChannels, length(Variable.nReflectors));
reSolution = cell(nChannels, length(Variable.nReflectors));

for iChannel = 1 : nChannels
    for iReflector = 1 : length(Variable.nReflectors)
        % * Get number of reflectors and define spatial correlation
        nReflectors = Variable.nReflectors(iReflector);
        corIrs = eye(nReflectors);

        % * Generate tap gains and delays
        [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
        [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
        [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

        % * Construct channels
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, subbandFrequency, fadingMode);

        % * Alternating optimization
        [reSample{iChannel, iReflector}, reSolution{iChannel, iReflector}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reSampleAvg = cell(1, length(Variable.nReflectors));
for iReflector = 1 : length(Variable.nReflectors)
    reSampleAvg{iReflector} = mean(cat(3, reSample{:, iReflector}), 3);
end
save('data/re_reflector.mat');

%% * R-E plots
figure('name', 'R-E region vs number of reflectors');
legendString = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    plot(reSampleAvg{iReflector}(1, :) / nSubbands, 1e6 * reSampleAvg{iReflector}(2, :));
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
