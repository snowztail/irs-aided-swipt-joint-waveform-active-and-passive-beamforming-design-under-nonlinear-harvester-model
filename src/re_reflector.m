clear; clc; setup; config_re_reflector;

%% ! R-E region vs number of IRS elements
reAoSample = cell(nChannels, length(Variable.nReflectors));
reLcSample = cell(nChannels, length(Variable.nReflectors));

reAoSolution = cell(nChannels, length(Variable.nReflectors));
reLcSolution = cell(nChannels, length(Variable.nReflectors));

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
        [directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
		[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

        % * Alternating optimization
        [reAoSample{iChannel, iReflector}, reAoSolution{iChannel, iReflector}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
        [reLcSample{iChannel, iReflector}, reLcSolution{iChannel, iReflector}] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    end
end

% * Average over channel realizations
reAoInstance = cell(1, length(Variable.nReflectors));
reLcInstance = cell(1, length(Variable.nReflectors));
flag = zeros(1, length(Variable.nReflectors));

for iReflector = 1 : length(Variable.nReflectors)
    reAoInstance{iReflector} = mean(cat(3, reAoSample{:, iReflector}), 3);
    reLcInstance{iReflector} = mean(cat(3, reLcSample{:, iReflector}), 3);
	flag(iReflector) = isempty(reAoInstance{iReflector}) || isempty(reLcInstance{iReflector});
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_reflector/re_reflector_%d.mat', iBatch));
end
