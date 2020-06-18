function [infoWaveform, powerWaveform, infoRatio, powerRatio, current, rate] = waveform_sdr(infoWaveform, powerWaveform, infoRatio, powerRatio, beta2, beta4, txPower, noisePower, rateConstraint, tolerance, compositeChannel, nCandidates, scaleFactor)
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
    %   - scaleFactor(\eta): scale factor on auxiliary variables to avoid precision problem
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



    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);
    % \boldsymbol{H}_{I/P}
    channelMatrix = compositeChannel * compositeChannel';
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
    % a^{(0)}, b^{(0)}, c^{(0)}
    infoProdVar = powerRatio * infoAuxiliary(nSubbands);
    powerProdVar = powerRatio * powerAuxiliary(nSubbands);
    powerAuxiliarySum = powerAuxiliary' * powerAuxiliary;

    % * SCA
    current_ = 0;
    isConverged = false;
    while ~isConverged
        % * Update solution, auxiliary, and SDR matrices
        infoRatio_ = infoRatio;
        powerRatio_ = powerRatio;
        infoMatrix_ = infoMatrix;
        infoAuxiliary_ = infoAuxiliary * scaleFactor;
        powerAuxiliary_ = powerAuxiliary * scaleFactor;
        infoProdVar_ = infoProdVar;
        powerProdVar_ = powerProdVar;
        powerAuxiliarySum_ = powerAuxiliarySum;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_precision high
            cvx_solver mosek
            variable infoMatrix(nSubbands, nSubbands) hermitian semidefinite;
            variable powerMatrix(nSubbands, nSubbands) hermitian semidefinite;
            variable infoRatio nonnegative;
            variable powerRatio nonnegative;
            variable infoProdVar nonnegative;
            variable powerProdVar nonnegative;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            expression signalPower(nSubbands, 1);
            expression sinr(nSubbands, 1);
            % t'_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
            end
            infoAuxiliary = infoAuxiliary * scaleFactor;
            powerAuxiliary = powerAuxiliary * scaleFactor;
            % \tilde{t}'_{I/P,0}
            infoProdLowerBound = (1 / 2) * (powerRatio + infoAuxiliary(nSubbands)) * (powerRatio_ + infoAuxiliary_(nSubbands)) ...
                - (1 / 4) * (powerRatio_ + infoAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (powerRatio - infoAuxiliary(nSubbands)) ^ 2;
            powerProdLowerBound = (1 / 2) * (powerRatio + powerAuxiliary(nSubbands)) * (powerRatio_ + powerAuxiliary_(nSubbands)) ...
                - (1 / 4) * (powerRatio_ + powerAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (powerRatio - powerAuxiliary(nSubbands)) ^ 2;
            % c
            powerAuxiliarySum = real(powerAuxiliary_' * powerAuxiliary);
            % g_p
            powerProdSumLowerBound = (1 / 2) * (powerAuxiliarySum + powerRatio) * (powerAuxiliarySum_ + powerRatio_) - (1 / 4) * (powerAuxiliarySum_ + powerRatio_) ^ 2 - (1 / 4) * (powerAuxiliarySum - powerRatio) ^ 2;
            % \tilde{z}'
            % currentLowerBound = (1 / 2) * beta2 * (infoProdLowerBound + powerProdLowerBound) ...
            %     + (3 / 8) * beta4 * ( ...
            %         8 * infoRatio_ * infoAuxiliary_(nSubbands) * infoProdLowerBound + 4 * powerRatio_ * powerProdSumLowerBound ...
            %         - (4 * powerRatio_ * infoAuxiliary_(nSubbands) ^ 2 + 2 * powerRatio_ * (powerAuxiliary_' * powerAuxiliary_)) * powerRatio ...
            %         - 4 * powerRatio_ ^ 2 * infoAuxiliary_(nSubbands) * infoAuxiliary(nSubbands) - 2 * real(powerRatio_ ^ 2 * powerAuxiliary_' * powerAuxiliary) ...
            %         + 2 * powerRatio_ ^ 2 * infoAuxiliary_(nSubbands) ^ 2 + powerRatio_ ^ 2 * (powerAuxiliary_' * powerAuxiliary_) ...
            %     ) ...
            %     + (3 / 2) * beta4 * ((1 / 2) * (infoProdVar + powerProdVar) * (infoProdVar_ + powerProdVar_) - (1 / 4) * (infoProdVar_ + powerProdVar_) ^ 2 - (1 / 4) * (infoProdVar - powerProdVar) ^ 2);
            currentLowerBound = (1 / 2) * beta2 * (infoProdLowerBound + powerProdLowerBound) ...
                + (3 / 8) * beta4 * ( ...
                    8 * infoRatio_ * infoAuxiliary_(nSubbands) * infoProdLowerBound + 4 * powerRatio_ * powerProdSumLowerBound ...
                    - (4 * powerRatio_ * infoAuxiliary_(nSubbands) ^ 2 + 2 * powerRatio_ * (powerAuxiliary_' * powerAuxiliary_)) * powerRatio ...
                    - 4 * powerRatio_ ^ 2 * infoAuxiliary_(nSubbands) * infoAuxiliary(nSubbands) - 2 * real(powerRatio_ ^ 2 * powerAuxiliary_' * powerAuxiliary) ...
                    + 2 * powerRatio_ ^ 2 * infoAuxiliary_(nSubbands) ^ 2 + powerRatio_ ^ 2 * (powerAuxiliary_' * powerAuxiliary_) ...
                );
            % g
            for iSubband = 1 : nSubbands
                signalPower(iSubband) = (1 / 2) * (infoRatio + infoMatrix(iSubband, iSubband)) * (infoRatio_ + infoMatrix_(iSubband, iSubband))...
                    - (1 / 4) * (infoRatio_ + infoMatrix_(iSubband, iSubband)) ^ 2 - (1 / 4) * (infoRatio - infoMatrix(iSubband, iSubband)) ^ 2;
            end
            % \gamma
            for iSubband = 1 : nSubbands
                sinr(iSubband) = signalPower(iSubband) * square_abs(compositeChannel(iSubband)) / noisePower;
            end
            maximize currentLowerBound;
            subject to
                (1 / 2) * (trace(infoMatrix) + trace(powerMatrix)) <= txPower;
                geo_mean(1 + sinr) >= 2 ^ (rateConstraint / nSubbands);
                powerRatio + infoRatio <= 1;
                % infoProdLowerBound >= infoProdVar;
                % (powerProdLowerBound) >= powerProdVar;
        cvx_end
        infoAuxiliary = infoAuxiliary / scaleFactor;
        powerAuxiliary = powerAuxiliary / scaleFactor;

        % * Update output current
        % z
        % current = ((1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
        %     + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
        %     + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));
        current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary));
        % R
        rate = 0;
        for iSubband = 1 : nSubbands
            rate = rate + log2(1 + infoRatio * infoMatrix(iSubband, iSubband) * square_abs(compositeChannel(iSubband)) / noisePower);
        end

        % * Test convergence
        (current - current_)
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end
    infoMatrix = full(infoMatrix);
    powerMatrix = full(powerMatrix);

    % * Recover rank-1 solution by randomization method
    [u1, sigma1] = eig(infoMatrix);
    [u2, sigma2] = eig(powerMatrix);
    current = 0;
    rate = 0;
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
        % R
        rate_ = 0;
        for iSubband = 1 : nSubbands
            rate_ = rate_ + log2(1 + infoRatio * infoMatrix(iSubband, iSubband) * square_abs(compositeChannel(iSubband)) / noisePower);
        end
        if current_ >= current && rate_ >= rateConstraint
            current = current_;
            rate = rate_;
            infoWaveform = infoWaveform_;
            powerWaveform = powerWaveform_;
        end
    end

end
