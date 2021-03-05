function [sample, solution] = re_sample_swipt_low_complexity(alpha, beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance)
    % Function:
    %   - sample R-E region by computing the output DC current and rate
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - nSamples (S): number of samples in R-E region
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - time sharing between WPT point achieved by SMF and WIT point achieved by WF
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    sample = zeros(2, nSamples);
    solution = cell(nSamples, 1);

	% * WIT point
    [sample(:, 1), solution{1}] = re_sample_wit_wf(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
	irs = solution{1}.irs;
    compositeChannel = solution{1}.compositeChannel;

	% * Remove rate constraint as R-E boundary is obtained by varying the waveform and splitting ratios
	rateConstraint = 0;

	for iSample = 2 : nSamples
		% * Update splitting ratio
		infoRatio = (nSamples - iSample) / (nSamples - 1);
		powerRatio = 1 - infoRatio;

		% TODO Update waveform ratio
		waveformRatio = powerRatio;

		isConverged = false;
		current_ = 0;
		eigRatio = [];
		while ~isConverged
			% * Design waveform by WF and SMF + MRT
			[~, infoAmplitude] = water_filling(compositeChannel, txPower, noisePower, waveformRatio);
			[~, ~, powerAmplitude] = scaled_matched_filter(alpha, beta2, beta4, compositeChannel, txPower, waveformRatio);
			[infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);

			% * Optimize IRS by SDR
			[irs, eigRatio(end + 1)] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates, tolerance);
			[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
			channelAmplitude = vecnorm(compositeChannel, 2, 2);

			% * Get R-E sample
			[rate] = water_filling(compositeChannel, txPower, noisePower, waveformRatio);
			[current] = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

			% * Check convergence
			isConverged = abs(current - current_) <= tolerance;
			current_ = current;
		end
		sample(:, iSample) = [rate; current];
        solution{iSample} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);
	end

end
