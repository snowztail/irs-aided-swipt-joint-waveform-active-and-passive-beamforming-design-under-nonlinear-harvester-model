function [infoAmplitude, powerAmplitude, current] = waveform_linear(beta2, beta4, channel, txPower)
    % Function:
    %   - optimize the information and power waveform amplitude to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
	%   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %
	% Output:
	%   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
	%   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
	%	- current (z): objective function to maximize output DC current
    %
	% Comment:
	%	- based on linear harvester model
    %   - for WPT, the optimal frequency allocation is ASS and the optimal spatial beamformer is MRT
    %   - solve SDR problem to obtain high-rank amplitude outer product matrices
    %   - use Gaussian randomization method to extract amplitude vectors
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

    % * Construct coefficient matrices
    % \boldsymbol{S}_{I/P}
    infoMatrix = vec(infoAmplitude) * vec(infoAmplitude)';
	powerMatrix = vec(powerAmplitude) * vec(powerAmplitude)';

    % \boldsymbol{A}_n
    channelMatrix = channelAmplitude * channelAmplitude';
	channelBlkDiag = block_diagonal(channelMatrix, 1, nSubbands);

    % * Initialize auxiliaries
    % t'_{I/P,n}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * powerMatrix);
    end

	% * Update output current
	current = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

end
