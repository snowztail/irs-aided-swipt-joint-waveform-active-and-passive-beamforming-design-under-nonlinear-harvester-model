function [current, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = waveform_ass(beta2, beta4, channel, txPower)
    % Function:
    %   - optimize the information and power waveform amplitude for WPT
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
	%   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
	% Output:
	%	- current (z): objective function to maximize output DC
	%   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
	%   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %
	% Comment:
	%	- based on linear harvester model
    %   - for linear WPT, the optimal frequency allocation is ASS and the optimal spatial beamformer is MRT
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Oct 20




    % * Get data
	[nSubbands] = size(channel, 1);

    % * Initialize algorithm
    % \boldsymbol{a}, \boldsymbol{s}_{I/P}
	channelAmplitude = vecnorm(channel, 2, 2);

	% * Adaptive single sine
	infoAmplitude = zeros(1, nSubbands) + eps;
	powerAmplitude = zeros(1, nSubbands) + eps;
	powerAmplitude(channelAmplitude == max(channelAmplitude)) = sqrt(2 * txPower);

	% * Splitting ratio
	infoRatio = eps;
	powerRatio = 1 - eps;

	% * Update output current
	[current] = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

end
