function [sample, solution] = re_sample_swipt_imperfect_csi(alpha, beta2, beta4, directChannel, cascadedChannel, imperfectCascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance)
    % Function:
    %   - sample R-E region by computing the output DC current and rate
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
	%   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
	%	- imperfectCascadedChannel (\hat{V}) [nReflectors * nTxs * nSubbands]: estimated AP-IRS-user concatenated channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - nSamples (S): number of samples in R-E region
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, perfect and imperfect composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - investigate the impact of imperfect CSIT on R-E behavior
    %
    % Author & Date: Yang (i@snowztail.com) - 12 Mar 21


	% * Solve problem based on imperfect CSIT
	[~, solution] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, imperfectCascadedChannel, txPower, noisePower, nCandidates, nSamples, tolerance);

	sample = zeros(2, nSamples);
	for iSample = 1 : nSamples
		% * Update composite channel
		solution{iSample}.imperfectCompositeChannel = solution{iSample}.compositeChannel;
		solution{iSample}.compositeChannel = composite_channel(directChannel, cascadedChannel, solution{iSample}.irs);

		% * Retrieve data
		channelAmplitude = vecnorm(solution{iSample}.compositeChannel, 2, 2);
		infoAmplitude = solution{iSample}.infoAmplitude;
		powerAmplitude = solution{iSample}.powerAmplitude;
		powerRatio = solution{iSample}.powerRatio;
		infoRatio = solution{iSample}.infoRatio;

		% * Calculate R-E pair
		sample(1, iSample) = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);
		sample(2, iSample) = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
	end

end
