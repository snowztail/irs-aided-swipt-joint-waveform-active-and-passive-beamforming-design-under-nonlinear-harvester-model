clear; clc; setup; config_reflector;

%% ! R-E region vs number of IRS elements
reSample = cell(nChannels, nCases);
reSolution = cell(nChannels, nCases);

parfor iChannel = 1 : nChannels
    for iReflector = 1 : nCases
        % * Get number of reflectors and define spatial correlation
        nReflectors = Variable.nReflectors(iReflector);
        corIrs = eye(nReflectors);

        % * Generate tap gains and delays
        [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
        [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
        [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

        % * Construct channels
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

        % * Alternating optimization
        [reSample{iChannel, iReflector}, reSolution{iChannel, iReflector}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reSampleAvg = cell(1, nCases);
for iReflector = 1 : nCases
    reSampleAvg{iReflector} = mean(cat(3, reSample{:, iReflector}), 3);
end
save('data/re_reflector.mat');

% %% * R-E plots
% figure('name', 'R-E region vs number of reflectors');
% legendString = cell(1, nCases);
% for iReflector = 1 : nCases
%     plot(reSampleAvg{iReflector}(1, :) / nSubbands, 1e6 * reSampleAvg{iReflector}(2, :));
%     legendString{iReflector} = sprintf('L = %d', Variable.nReflectors(iReflector));
%     hold on;
% end
% hold off;
% grid minor;
% legend(legendString);
% xlabel('Per-subband rate [bps/Hz]');
% ylabel('Average output DC current [\muA]');
% ylim([0 inf]);
% savefig('plots/re_reflector.fig');
