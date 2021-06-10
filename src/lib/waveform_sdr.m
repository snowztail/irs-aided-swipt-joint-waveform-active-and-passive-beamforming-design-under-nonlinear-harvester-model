function [current, infoAmplitude, powerAmplitude] = waveform_sdr(beta2, beta4, channel, infoAmplitude, powerAmplitude, txPower, nCandidates, tolerance)
    % Function:
    %   - optimize the information and power waveform amplitude to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain (in the previous iteration)
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain (in the previous iteration)
    %   - txPower (P): average transmit power budget
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
	%	- current (z): objective function to maximize output DC current
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
	%   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %
    % Comment:
    %   - suitable for WPT, the result inlines with GP but with much lower complexity
    %   - solve SDR problem to obtain high-rank amplitude outer product matrices
    %   - use Gaussian randomization method to extract amplitude vectors
    %
    % Author & Date: Yang (i@snowztail.com) - 31 May 20


    % * Get data
    [nSubbands] = size(channel, 1);

    % * Initialize algorithm
    % \boldsymbol{a}, \boldsymbol{s}_{I/P}
    channelAmplitude = vecnorm(channel, 2, 2);

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

    % * SCA
    current_ = 0;
    isConverged = false;
    while ~isConverged
        % * Update auxiliaries
        infoAuxiliary_ = infoAuxiliary;
        powerAuxiliary_ = powerAuxiliary;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            variable infoMatrix(nSubbands, nSubbands) hermitian semidefinite;
            variable powerMatrix(nSubbands, nSubbands) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            % t'_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * infoMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * powerMatrix);
            end
            % \tilde{z}
            currentSca = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + 2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_) ...
				+ (3 / 2) * beta4 * (infoAuxiliary(nSubbands) * powerAuxiliary_(nSubbands) + powerAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) * powerAuxiliary_(nSubbands));
            maximize currentSca;
            subject to
                (1 / 2) * (trace(infoMatrix) + trace(powerMatrix)) <= txPower;
        cvx_end

        % * Update output current
        current = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Test convergence
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
    end
    infoMatrix = full(infoMatrix);
    powerMatrix = full(powerMatrix);

    % * Recover rank-1 solution by randomization method
    [infoU, infoSigma] = eig(infoMatrix);
    [powerU, powerSigma] = eig(powerMatrix);
    current = 0;
	for iCandidate = 1 : nCandidates
        infoAmplitudeCandidate = infoU * infoSigma ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        powerAmplitudeCandidate = powerU * powerSigma ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        infoMatrix = infoAmplitudeCandidate * infoAmplitudeCandidate';
        powerMatrix = powerAmplitudeCandidate * powerAmplitudeCandidate';

        % * Update auxiliaries
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * infoMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * powerMatrix);
        end

        % * Compute current
        currentCandidate = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Choose best candidate
        if currentCandidate >= current
            current = currentCandidate;
            infoAmplitude = transpose(infoAmplitudeCandidate);
            powerAmplitude = transpose(powerAmplitudeCandidate);
        end
	end

end
