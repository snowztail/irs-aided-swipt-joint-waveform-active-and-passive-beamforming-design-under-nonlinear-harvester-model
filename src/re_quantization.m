clear; clc; setup; config_re_quantization;

%% ! R-E region vs quantization bits
reNoIrsSample = cell(nChannels, 1);
reIrsSample = cell(nChannels, 1);
reQuantizedSample = cell(nChannels, length(Variable.nQuantizeBits));

reNoIrsSolution = cell(nChannels, 1);
reIrsSolution = cell(nChannels, 1);
reQuantizedSolution = cell(nChannels, length(Variable.nQuantizeBits));

for iChannel = 1 : nChannels
	% * Generate tap gains and delays
	[directTapGain, directTapDelay] = tap_tgn(corTx, corRx, 'nlos');
	[incidentTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, 'nlos');
	[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, 'nlos');

	% * Construct channels
	[directChannel] = channel_response(directTapGain, directTapDelay, directDistance, rxGain, subbandFrequency, fadingMode);
	[incidentChannel] = channel_response(incidentTapGain, incidentTapDelay, incidentDistance, irsGain, subbandFrequency, fadingMode);
	[reflectiveChannel] = channel_response(reflectiveTapGain, reflectiveTapDelay, reflectiveDistance, rxGain, subbandFrequency, fadingMode);
	[cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel);

	% * R-E region without IRS
	[reNoIrsSample{iChannel}, reNoIrsSolution{iChannel}] = re_sample_swipt_benchmark(alpha, beta2, beta4, directChannel, txPower, noisePower, nSamples, tolerance);

	% * R-E region with IRS
	[reIrsSample{iChannel}, reIrsSolution{iChannel}] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

    for iBit = 1 : length(Variable.nQuantizeBits)
		% * Update number of quantization bits
        nQuantizeBits = Variable.nQuantizeBits(iBit);

        % * R-E region with quantized IRS
        [reQuantizedSample{iChannel, iBit}, reQuantizedSolution{iChannel, iBit}] = re_sample_swipt_quantized_irs(beta2, beta4, directChannel, cascadedChannel, noisePower, nQuantizeBits, reIrsSolution{iChannel});
    end
end

% * Average over channel realizations
reQuantizedInstance = cell(1, length(Variable.nQuantizeBits));
flag = zeros(1, length(Variable.nQuantizeBits));

reNoIrsInstance = mean(cat(3, reNoIrsSample{:}), 3);
reIrsInstance = mean(cat(3, reIrsSample{:}), 3);
for iBit = 1 : length(Variable.nQuantizeBits)
    reQuantizedInstance{iBit} = mean(cat(3, reQuantizedSample{:, iBit}), 3);
	flag(iBit) = isempty(reNoIrsInstance) || isempty(reIrsInstance) || isempty(reQuantizedInstance{iBit});
end

% * Save batch data
if ~sum(flag(:))
	save(sprintf('data/re_quantization/re_quantization_%d.mat', iBatch));
end
