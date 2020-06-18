function [irs, current, rate] = irs_flat(irs, beta2, beta4, noisePower, rateConstraint, tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio, nCandidates)
    % Function:
    %   - optimize the IRS reflection coefficients to maximize the R-E region
    %   - compute the output DC current and user rate
    %
    % Input:
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients (in the previous iteration)
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - concatVector (M) [(nReflectors + 1) * nSubbands]: concatenated channel vector
    %   - concatMatrix (R_n) [(nReflectors + 1) * (nReflectors + 1)]: rate SDR matrix
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (in the previous iteration)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (in the previous iteration)
    %   - infoRatio (\bar{\rho}): information splitting ratio (in the previous iteration)
    %   - powerRatio (\rho): power splitting ratio (in the previous iteration)
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %
    % Output:
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - current: maximum achievable output DC current
    %   - rate: maximum achievable user rate
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);
    nReflectors = size(concatMatrix{1}, 1) - 1;
    % \bar{\boldsymbol{\phi}}^{(0)}
    irs = [irs; 1];
    % \bar{\boldsymbol{\Phi}}^{(0)}
    irsMatrix = irs * irs';
    % \boldsymbol{W}_{I/P}
    infoMatrix = infoWaveform * infoWaveform';
    powerMatrix = powerWaveform * powerWaveform';
    % \boldsymbol{C}_{I/P,n}
    infoCoefMatrix = cell(2 * nSubbands - 1, 1);
    powerCoefMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoCoefMatrix{iSubband + nSubbands} = concatVector * diag(diag(conj(infoMatrix), iSubband), iSubband) * concatVector';
        powerCoefMatrix{iSubband + nSubbands} = concatVector * diag(diag(conj(powerMatrix), iSubband), iSubband) * concatVector';
    end
    infoCoefMatrix{nSubbands} = hermitianize(infoCoefMatrix{nSubbands});
    powerCoefMatrix{nSubbands} = hermitianize(powerCoefMatrix{nSubbands});
    % t_{I/P,n}^{(0)}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
    end
    infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
    powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

    % * SCA
    current_ = 0;
    isConverged = false;
    while ~isConverged
        % * Update auxiliary and SDR matrices
        % t_{I/P,n}^{(i-1)}
        infoAuxiliary_ = infoAuxiliary;
        powerAuxiliary_ = powerAuxiliary;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            expression sinr(nSubbands, 1);
            % t_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
            end
            infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
            powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

            % \tilde{z}
            currentLowerBound = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + (2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_)) ...
                + (3 / 2) * beta4 * powerRatio ^ 2 * ((1 / 2) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) - (1 / 4) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (infoAuxiliary(nSubbands) - powerAuxiliary(nSubbands)) ^ 2);
            % \gamma
            for iSubband = 1 : nSubbands
                sinr(iSubband) = infoRatio * square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower;
            end
            maximize currentLowerBound;
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
                geo_mean(1 + sinr) >= 2 ^ (rateConstraint / nSubbands);
        cvx_end

        % * Update output current
        % z
        current = real((1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end
    irsMatrix = full(irsMatrix);

    % * Recover rank-1 solution by randomization method
    [u, sigma] = eig(irsMatrix);
    current = 0;
    rate = 0;
    for iCandidate = 1 : nCandidates
        % \bar{\boldsymbol{\phi}}_q
        irs_ = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
        % \bar{\boldsymbol{\Phi}}_q
        irsMatrix = irs_ * irs_';
        % t_{I/P,n}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
        end
        % z
        current_ = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);
        % R
        rate_ = 0;
        for iSubband = 1 : nSubbands
            rate_ = rate_ + log2(1 + infoRatio * square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower);
        end
        if current_ >= current && rate_ >= rateConstraint
            current = current_;
            rate = rate_;
            irs = irs_;
        end
    end
    % \boldsymbol{\phi}
    irs = irs(1 : nReflectors) / irs(end);

end
