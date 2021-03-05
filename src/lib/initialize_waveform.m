function [infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize_waveform(alpha, beta2, beta4, channel, txPower, noisePower)
    % Function:
    %   - initialize waveform and splitting ratio for GP algorithm
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - the iterative GP algorithm is sensitive to initialization (only converges to stationary points)
    %   - GP requires nonzero entries thus use eps to replace zero
    %   - cvx crashes when initial rate is far from constraint
    %       - use scaled matched filter to initialize power waveform with power P
    %       - use water-filling + MRT to initialize information waveform with power P
    %       - initialize both splitting ratio to 1
	%		- regulated over iterations
	%	- if the number of subbands is too small, it is unnecessary to use multisine waveform
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Jul 20


	% * Get data
	nSubbands = size(channel, 1);

	% * Initialize modulated waveform by WF
	[~, infoAmplitude] = water_filling(channel, txPower, noisePower);

	% * If necessary, initialize multisine waveform by SMF
	if nSubbands <= 2
		powerAmplitude = zeros(1, nSubbands) + eps;
	else
		[~, ~, powerAmplitude] = scaled_matched_filter(alpha, beta2, beta4, channel, txPower);
	end

	infoRatio = 1 - eps;
    powerRatio = 1 - eps;

end
