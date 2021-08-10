function [current, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = scaled_matched_filter(alpha, beta2, beta4, channel, txPower, waveformRatio)
    % Function:
    %   - power allocation based on scaled matched filter for OFDM channels
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
	%	- waveformRatio (\delta): balancing ratio for modulated and multisine waveforms
    %
    % Output:
	%	- current (z): objective function to maximize output DC
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - for MISO OFDM channels
	%	- convert waveform amplitude to peak value by multiplying sqrt(2)
	%	- assign no power to modulated waveform and set splitting ratio to 1
    %
    % Author & Date: Yang (i@snowztail.com) - 05 Mar 21


	if nargin == 5
		waveformRatio = 1;
	end
    % * Get data
    nSubbands = size(channel, 1);
    channelAmplitude = vecnorm(channel, 2, 2);

    % * SMF power allocation based on equivalent channel strength
	powerAmplitude = sqrt(2 * waveformRatio * txPower / sum(channelAmplitude .^ (2 * alpha))) * channelAmplitude' .^ alpha;

	% * Assign modulated waveform and splitting ratio
    infoAmplitude = zeros(1, nSubbands) + eps;
    infoRatio = eps;
    powerRatio = 1 - eps;

	% * Update output current
	[current] = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

end
