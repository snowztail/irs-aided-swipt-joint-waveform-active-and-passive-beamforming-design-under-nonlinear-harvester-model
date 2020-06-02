function [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint, tolerance, compositeChannel, nCandidates)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %   - compute the output DC current and user rate
    %
    % Input:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (in the previous iteration)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (in the previous iteration)
    %   - infoRatio (\bar{\rho}): information splitting ratio (in the previous iteration)
    %   - powerRatio (\rho): power splitting ratio (in the previous iteration)
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - txPower (P): transmit power constraint
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: the composite channel
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - current: maximum achievable output DC current
    %   - rate: maximum achievable user rate
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank waveform outer product matrices
    %   - use Gaussian randomization method to extract waveform vectors
    %
    % Author & Date: Yang (i@snowztail.com) - 31 May 20



    % * Construct current SDR matrices
    nSubbands = size(infoWaveform, 1);
    % \boldsymbol{H}_{I/P}
    channelMatrix = compositeChannel * compositeChannel';
    % \boldsymbol{H}_{I/P,n}
    channelCoefMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        channelCoefMatrix{iSubband + nSubbands} = diag(diag(channelMatrix, iSubband), iSubband);
    end

    % * SCA
    current_ = 0;
    isConverged = false;
    infoMatrix = infoWaveform * infoWaveform';
    powerMatrix = powerWaveform * powerWaveform';
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    % t'_{I/P,n}
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
    end

    while ~isConverged
        % \boldsymbol{A}^{(i)}
        % infoCoefMatrix = (1 / 2 * beta2) * conj(channelCoefMatrix{nSubbands}) + (3 / 2 * beta4) * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) * conj(channelCoefMatrix{nSubbands});
        infoCoefMatrix = (1 / 2 * beta2) * conj(channelCoefMatrix{nSubbands}) + (3 / 2 * beta4) * (infoAuxiliary(nSubbands)) * conj(channelCoefMatrix{nSubbands});
        powerCoefMatrix = (1 / 2 * beta2) * conj(channelCoefMatrix{nSubbands});
        for iSubband = - nSubbands + 1 : nSubbands - 1
            powerCoefMatrix = powerCoefMatrix + (3 / 8 * beta4) ...
                * (conj(powerAuxiliary(iSubband + nSubbands)) * conj(channelCoefMatrix{iSubband + nSubbands}) + powerAuxiliary(iSubband + nSubbands) * transpose(channelCoefMatrix{iSubband + nSubbands}));
        end
        % powerCoefMatrix = powerCoefMatrix + (3 / 2 * beta4) * infoAuxiliary(nSubbands) * conj(channelCoefMatrix{nSubbands});

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            variable infoMatrix(nSubbands, nSubbands) hermitian semidefinite;
            variable powerMatrix(nSubbands, nSubbands) hermitian semidefinite;
            % variable powerRatio;
            expression current;
            expression rate;
            % \tilde(z)
            current = powerRatio * (trace(infoCoefMatrix * infoMatrix) + trace(powerCoefMatrix * powerMatrix));
            % current = powerRatio * (trace(infoCoefMatrix * infoMatrix) + trace(powerCoefMatrix * powerMatrix)) ...
            %     - (3 / 8 * beta4) * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            %     - (3 / 2 * beta4) * (infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));
            % R
            for iSubband = 1 : nSubbands
                rate = rate + log(1 + (1 - powerRatio) * infoMatrix(iSubband, iSubband) * square_abs(compositeChannel(iSubband)) / noisePower) / log(2);
            end
            maximize current;
            subject to
                (1 / 2) * (trace(infoMatrix) + trace(powerMatrix)) <= txPower;
                rate >= rateConstraint;
        cvx_end
        (1 / 2) * (trace(infoMatrix) + trace(powerMatrix))

        % * Update auxiliary variables and output current
        % t'_{I/P,n}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
        end
        % z
        current = (1 / 2 * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + 3 / 8 * beta4 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            * 3 / 2 * beta4 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance || current == 0 || isnan(current);
        current_ = current;
    end
    infoMatrix = full(infoMatrix);
    powerMatrix = full(powerMatrix);

    % * Recover rank-1 solution by randomization method
    [u1, sigma1] = eig(infoMatrix);
    [u2, sigma2] = eig(powerMatrix);
    current = 0;
    for iCandidate = 1 : nCandidates
        infoWaveform_ = u1 * sigma1 ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        powerWaveform_ = u2 * sigma2 ^ (1 / 2) * exp(1i * 2 * pi * rand(nSubbands, 1));
        infoMatrix = infoWaveform_ * infoWaveform_';
        powerMatrix = powerWaveform_ * powerWaveform_';
        % t'_{I/P,n}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerRatio * conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
        end
        % z
        current_ = (1 / 2 * beta2 * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + 3 / 8 * beta4 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            * 3 / 2 * beta4 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));
        % R
        rate_ = 0;
        for iSubband = 1 : nSubbands
            rate_ = rate_ + log2(1 + infoRatio * square_abs(infoWaveform_(iSubband)) * square_abs(compositeChannel(iSubband)) / noisePower);
        end
        if current_ > current && rate_ >= rateConstraint
            current = current_;
            rate = rate_;
            infoWaveform = infoWaveform_;
            powerWaveform = powerWaveform_;
        end
    end

end
