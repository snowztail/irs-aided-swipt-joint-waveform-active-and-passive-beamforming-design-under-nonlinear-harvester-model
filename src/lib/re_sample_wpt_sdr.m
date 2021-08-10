function [sample, solution] = re_sample_wpt_sdr(alpha, beta2, beta4, directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize average output DC
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
    %   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
	%   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - solve both problems by SDR
	%	- for waveform design, the result inlines with GP but with much lower complexity
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    % * Get data
	[nReflectors, ~, ~] = size(cascadedChannel);

    % * Initialize IRS and composite channel
    irs = exp(1i * 2 * pi * rand(nReflectors, 1));
    [compositeChannel] = composite_channel(directChannel, cascadedChannel, irs);

    % * Initialize waveform and splitting ratio
	[~, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = scaled_matched_filter(alpha, beta2, beta4, compositeChannel, txPower);
	[infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);

    % * AO
    isConverged = false;
    current_ = 0;
	rateConstraint = 0;
	eigRatio = [];
    while ~isConverged
		[irs, eigRatio(end + 1)] = irs_sdr(beta2, beta4, directChannel, cascadedChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates, tolerance);
		[compositeChannel] = composite_channel(directChannel, cascadedChannel, irs);
		[current, infoAmplitude, powerAmplitude] = waveform_sdr(beta2, beta4, compositeChannel, infoAmplitude, powerAmplitude, txPower, nCandidates, tolerance);
		[infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
    end

	sample = [eps; current];
	solution = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);

end
