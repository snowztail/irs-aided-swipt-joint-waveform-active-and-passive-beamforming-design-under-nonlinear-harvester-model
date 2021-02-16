clear; clc; setup; config_scaling_tx;

%% ! WIT/WPT vs number of IRS elements
rateSample = zeros(nChannels, length(Variable.nTxs));
witSolution = cell(nChannels, length(Variable.nTxs));

currentLinearSample = zeros(nChannels, length(Variable.nTxs));
currentNonlinearSample = zeros(nChannels, length(Variable.nTxs));
wptLinearSolution = cell(nChannels, length(Variable.nTxs));
wptNonlinearSolution = cell(nChannels, length(Variable.nTxs));

for iChannel = 1 : nChannels
    for iTx = 1 : length(Variable.nTxs)
        % * Get number of transmit antennas and define spatial correlation
        nTxs = Variable.nTxs(iTx);
        corTx = eye(nTxs);

        % * Generate tap gains and delays
        [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
        [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
        [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

        % * Construct channels
        [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
        [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
        [reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

		% * WIT
		[rateSample(iChannel, iTx), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		witSolution{iChannel, iTx} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio eigRatio;

		% * WPT by linear harvester model
		[currentLinearSample(iChannel, iTx), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio] = wpt_linear(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		wptLinearSolution{iChannel, iTx} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio eigRatio;

		% * WPT by nonlinear harvester model
		[currentNonlinearSample(iChannel, iTx), irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio] = wpt(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		wptNonlinearSolution{iChannel, iTx} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);
		clearvars irs compositeChannel infoAmplitude powerAmplitude infoRatio powerRatio eigRatio;
    end
end

% * Average over channel realizations
rateInstance = mean(rateSample, 1);
currentLinearInstance = mean(currentLinearSample, 1);
currentNonlinearInstance = mean(currentNonlinearSample, 1);

% * Save batch data
save(sprintf('data/scaling_tx/scaling_tx_%d.mat', iBatch));
