function [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint, tolerance, infoRatio, powerRatio, noisePower, channel, infoWaveform, powerWaveform)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %   - based on SDR
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - txPower (P): transmit power constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (in the previous iteration)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (in the previous iteration)
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank waveform outer product matrices
    %   - use Gaussian randomization method to extract waveform vectors
    %
    % Author & Date: Yang (i@snowztail.com) - 31 May 20



    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);
    % \boldsymbol{H}_{I/P}
    channelMatrix = channel * channel';
    % \boldsymbol{H}_{I/P,n}
    channelCoefMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        channelCoefMatrix{iSubband + nSubbands} = diag(diag(channelMatrix, iSubband), iSubband);
    end
    % \boldsymbol{W}_{I/P}^{(0)}
    infoMatrix = infoWaveform * infoWaveform';
    powerMatrix = powerWaveform * powerWaveform';
    % t'_{I/P,n}^{(0)}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
    end

    % * SCA
    current_ = 0;
    isConverged = false;
    while ~isConverged
        % * Update auxiliary variables
        % t_{I/P,n}^{(i-1)}
        infoAuxiliary_ = infoAuxiliary;
        powerAuxiliary_ = powerAuxiliary;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            cvx_precision high
            variable infoMatrix(nSubbands, nSubbands) hermitian semidefinite;
            variable powerMatrix(nSubbands, nSubbands) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            expression snr(nSubbands, 1);
            % t'_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
            end
            % \tilde{z}'
            currentLb = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + 2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_) ...
                + (3 / 2) * beta4 * powerRatio ^ 2 * ((1 / 2) * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) - (1 / 4) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (infoAuxiliary(nSubbands) - powerAuxiliary(nSubbands)) ^ 2);
            % \gamma
            for iSubband = 1 : nSubbands
                snr(iSubband) = infoRatio * infoMatrix(iSubband, iSubband) * abs(channel(iSubband)) ^ 2 / noisePower;
            end
            maximize currentLb;
            subject to
                (1 / 2) * (trace(infoMatrix) + trace(powerMatrix)) <= txPower;
                geo_mean(1 + snr) >= 2 ^ (rateConstraint / nSubbands);
        cvx_end

        % * Update output current
        % z
        current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance || current == 0;
        current_ = current;
    end
    infoMatrix = full(infoMatrix);
    powerMatrix = full(powerMatrix);

    % * Recover rank-1 solution by randomization method
    [u1, sigma1] = eig(infoMatrix);
    [u2, sigma2] = eig(powerMatrix);
    current = 0;
    for iCandidate = 1 : nCandidates
        % \boldsymbol{w}_{I/P,q}
        infoWaveform_ = u1 * sigma1 ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        powerWaveform_ = u2 * sigma2 ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        % \boldsymbol{W}_{I/P,q}
        infoMatrix = infoWaveform_ * infoWaveform_';
        powerMatrix = powerWaveform_ * powerWaveform_';
        % t'_{I/P,n}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
        end
        % z
        current_ = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);
        % \boldsymbol{w}_{I/P}^\star
        if current_ >= current
            current = current_;
            infoWaveform = infoWaveform_;
            powerWaveform = powerWaveform_;
        end
    end

end
