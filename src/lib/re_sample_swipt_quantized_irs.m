function [sample, solution] = re_sample_swipt_quantized_irs(beta2, beta4, directChannel, cascadedChannel, noisePower, nQuantizeBits, solution)
    % Function:
    %   - sample R-E region by computing the output DC current and rate
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
	%   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %   - noisePower (\sigma_n^2): average noise power
    %	- nQuantizeBits (b): number of bits used for quantization
    %   - solution: unquantized IRS reflection coefficient, perfect and imperfect composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: quantized and unquantized IRS reflection coefficient, composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - investigate the impact of IRS quantization on R-E behavior
    %
    % Author & Date: Yang (i@snowztail.com) - 12 Mar 21


	[nSamples, ~] = size(solution);
	sample = zeros(2, nSamples);

	for iSample = 1 : nSamples
		% * Quantize IRS
		[solution{iSample}.quantizedIrs] = quantized_irs(solution{iSample}.irs, nQuantizeBits);

		% * Update composite channel
		solution{iSample}.quantizedCompositeChannel = composite_channel(directChannel, cascadedChannel, solution{iSample}.quantizedIrs);

		% * Retrieve data
		channelAmplitude = vecnorm(solution{iSample}.quantizedCompositeChannel, 2, 2);
		infoAmplitude = solution{iSample}.infoAmplitude;
		powerAmplitude = solution{iSample}.powerAmplitude;
		powerRatio = solution{iSample}.powerRatio;
		infoRatio = solution{iSample}.infoRatio;

		% * Calculate R-E pair
		sample(1, iSample) = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);
		sample(2, iSample) = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
	end

end
