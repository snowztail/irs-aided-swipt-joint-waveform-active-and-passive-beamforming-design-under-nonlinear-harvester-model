function [infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize_waveform(channel, txPower, noisePower)
    % Function:
    %   - initialize waveform and splitting ratio for GP algorithm
    %
    % Input:
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
    %   - the iterative GP algorithm is sensitive to initialization (only converges to stationary point)
    %   - GP requires nonzero entries thus use eps to replace zero
    %   - cvx crashes when initial rate is far from constraint
    %       - use scaled matched filter to initialize power waveform with power P
    %       - use water-filling + MRT to initialize information waveform with power P
    %       - initialize both splitting ratio to 1
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Jul 20


    % * Get equivalent channel gain
    channelAmplitude = vecnorm(channel, 2, 2);

    % * Initialize algorithm
    infoRatio = 1;
    powerRatio = 1;
	[~, infoAmplitude] = water_filling(channel, txPower, noisePower);

	% * If number of subbands is too small, it is unnecessary to use power waveform
	nSubbands = size(channel, 1);
	if nSubbands <= 2
		powerAmplitude = zeros(1, nSubbands) + eps;
	else
		powerAmplitude = sqrt(2 * txPower) * channelAmplitude' ./ norm(channelAmplitude);
	end

end
