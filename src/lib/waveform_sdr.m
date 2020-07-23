function [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, channel, infoWaveform, powerWaveform, txPower, nCandidates, tolerance)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform (in the previous iteration)
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform (in the previous iteration)
    %   - txPower (P): transmit power constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %
    % Comment:
    %   - suitable for WPT
    %   - solve SDR problem to obtain high-rank waveform outer product matrices
    %   - use Gaussian randomization method to extract waveform vectors
    %
    % Author & Date: Yang (i@snowztail.com) - 31 May 20


    % * Get data
    [nSubbands, nTxs, ~] = size(channel);

    % * Construct coefficient matrices
    % \boldsymbol{W}_{I/P}, \boldsymbol{W}_{I/P,n}
    infoMatrix = vec(infoWaveform) * vec(infoWaveform)';
    powerMatrix = vec(powerWaveform) * vec(powerWaveform)';

    % \boldsymbol{H}_n
    channelMatrix = vec(channel') * vec(channel')';
    channelBlkDiag = block_diagonal(channelMatrix, nTxs, nSubbands);

    % * Initialize auxiliaries
    % t'_{I/P,n}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * powerMatrix);
    end
    infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
    powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

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
            % cvx_precision high
            variable infoMatrix(nTxs * nSubbands, nTxs * nSubbands) hermitian semidefinite;
            variable powerMatrix(nTxs * nSubbands, nTxs * nSubbands) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            % t'_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * infoMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * powerMatrix);
            end
            infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
            powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));
            % \tilde{z}
            currentSca = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 ^ 2 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + 2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_) ...
                + (3 / 2) * beta4 ^ 2 * ((1 / 2) * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) - (1 / 4) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (infoAuxiliary(nSubbands) - powerAuxiliary(nSubbands)) ^ 2);
            maximize currentSca;
            subject to
                (1 / 2) * (trace(infoMatrix) + trace(powerMatrix)) <= txPower;
        cvx_end

        % * Update output current
        current = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end
    infoMatrix = full(infoMatrix);
    powerMatrix = full(powerMatrix);

    % * Recover rank-1 solution by randomization method
    [u1, sigma1] = eig(infoMatrix);
    [u2, sigma2] = eig(powerMatrix);
    current = 0;
    for iCandidate = 1 : nCandidates
        infoWaveformCandidate = u1 * sigma1 ^ (1 / 2) * exp(1i * 2 * pi * rand(nTxs * nSubbands, 1));
        powerWaveformCandidate = u2 * sigma2 ^ (1 / 2) * exp(1i * 2 * pi * rand(nTxs * nSubbands, 1));
        infoMatrix = infoWaveformCandidate * infoWaveformCandidate';
        powerMatrix = powerWaveformCandidate * powerWaveformCandidate';

        % * Update auxiliaries
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * infoMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(conj(channelBlkDiag{iSubband + nSubbands}) * powerMatrix);
        end
        infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
        powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

        % * Compute current
        currentCandidate = (1 / 2) * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Choose best candidate
        if currentCandidate >= current
            current = currentCandidate;
            infoWaveform = reshape(infoWaveformCandidate, [nTxs, nSubbands]);
            powerWaveform = reshape(powerWaveformCandidate, [nTxs, nSubbands]);
        end
    end

end
