function [current] = current_sdr(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio)
    % Function:
    %   - formulate output DC current as a function of waveform amplitudes
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channelAmplitude (a) [nSubbands * 1]: channel amplitude response
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude on information carriers
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude on power carriers
    %   - powerRatio (\rho): power splitting ratio
    %
    % Output:
    %	- current (z): objective function to maximize output DC current
    %
    % Comment:
    %   - based on SDR expression
    %
    % Author & Date: Yang (i@snowztail.com) - 05 Mar 21


	% * Get data
	nSubbands = size(channelAmplitude, 1);

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
	current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
		+ (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
		+ (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

end
