clear; clc; setup; config_irs;

%% ! R-E region for fixed and adaptive IRS
reAdaptiveSample = cell(nChannels, 1);
reIdealSample = cell(nChannels, 1);
reWitSample = cell(nChannels, 1);
reWptSample = cell(nChannels, 1);
reNoIrsSample = cell(nChannels, 1);

reAdaptiveSolution = cell(nChannels, 1);
reIdealSolution = cell(nChannels, 1);
reWitSolution = cell(nChannels, 1);
reWptSolution = cell(nChannels, 1);
reNoIrsSolution = cell(nChannels, 1);

for iChannel = 1 : nChannels
    % * Generate tap gains and delays
    [directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
    [incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
    [reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

    % * Construct channels
    [directChannel] = frequency_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
    [incidentChannel] = frequency_response(incidentTapGain, incidentTapDelay, incidentDistance, rxGain, subbandFrequency, fadingMode);
	[reflectiveChannel] = frequency_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);

    % * Adaptive IRS and waveform design
    [reAdaptiveSample{iChannel}, reAdaptiveSolution{iChannel}] = re_sample(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance);
    compositeChannelWit = reAdaptiveSolution{iChannel}{1}.compositeChannel;
    compositeChannelWpt = reAdaptiveSolution{iChannel}{end}.compositeChannel;

	% * Upper bound by frequency-selective IRS
	[compositeChannelIdeal] = composite_channel_ideal(directChannel, incidentChannel, reflectiveChannel);
	[reIdealSample{iChannel}, reIdealSolution{iChannel}] = re_sample_reference(beta2, beta4, compositeChannelIdeal, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization with fixed WIT-optimized IRS
    [reWitSample{iChannel}, reWitSolution{iChannel}] = re_sample_reference(beta2, beta4, compositeChannelWit, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization with fixed WPT-optimized IRS
    [reWptSample{iChannel}, reWptSolution{iChannel}] = re_sample_reference(beta2, beta4, compositeChannelWpt, txPower, noisePower, nSamples, tolerance);

    % * Waveform optimization without IRS
    [reNoIrsSample{iChannel}, reNoIrsSolution{iChannel}] = re_sample_reference(beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);
end

% * Average over channel realizations
reInstance{1} = mean(cat(3, reAdaptiveSample{:}), 3);
reInstance{2} = mean(cat(3, reIdealSample{:}), 3);
reInstance{3} = mean(cat(3, reWitSample{:}), 3);
reInstance{4} = mean(cat(3, reWptSample{:}), 3);
reInstance{5} = mean(cat(3, reNoIrsSample{:}), 3);

% * Save batch data
save(sprintf('data/re_irs_%d.mat', iBatch));
