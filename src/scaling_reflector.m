clear; clc; setup; config_scaling_reflector;

%% ! WIT/WPT vs number of IRS elements
rateSample = zeros(nChannels, length(Variable.nReflectors));
witSolution = cell(nChannels, length(Variable.nReflectors));

currentLinearSample = zeros(nChannels, length(Variable.nReflectors));
currentNonlinearSample = zeros(nChannels, length(Variable.nReflectors));
wptLinearSolution = cell(nChannels, length(Variable.nReflectors));
wptNonlinearSolution = cell(nChannels, length(Variable.nReflectors));

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
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

		% * WIT
		[rateSample(iChannel, iReflector), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		witSolution{iChannel, iReflector} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio;

		% * WPT by linear harvester model
		[currentLinearSample(iChannel, iReflector), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = wpt_linear(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		wptLinearSolution{iChannel, iReflector} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio;

		% * WPT by nonlinear harvester model
		[currentNonlinearSample(iChannel, iReflector), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = wpt(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		wptNonlinearSolution{iChannel, iReflector} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio;
    end
end

% * Average over channel realizations
rateInstance = mean(rateSample, 1);
currentLinearInstance = mean(currentLinearSample, 1);
currentNonlinearInstance = mean(currentNonlinearSample, 1);

% * Save batch data
save(sprintf('data/scaling_reflector/scaling_reflector_%d.mat', iBatch));
