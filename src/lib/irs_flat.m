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



    % * Construct current SDR matrices
    nSubbands = size(infoWaveform, 1);
    nReflectors = size(concatMatrix{1}, 1) - 1;
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

    % * SCA
    current_ = 0;
    isConverged = false;
    irs = [irs; 1];
    irsMatrix = irs * irs';
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    while ~isConverged
        % t_{I/P,n}^{(i-1)}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
        end
        % \boldsymbol{A}^{(i)}
        coefMatrix = (1 / 2 * beta2) * powerRatio * (infoCoefMatrix{nSubbands} + powerCoefMatrix{nSubbands});
        for iSubband = - nSubbands + 1 : nSubbands - 1
            coefMatrix = coefMatrix + (3 / 8 * beta4) * powerRatio ^ 2 ...
                * (2 * (conj(infoAuxiliary(iSubband + nSubbands)) * infoCoefMatrix{iSubband + nSubbands} + infoAuxiliary(iSubband + nSubbands) * ctranspose(infoCoefMatrix{iSubband + nSubbands})) ...
                + (conj(powerAuxiliary(iSubband + nSubbands)) * powerCoefMatrix{iSubband + nSubbands} + powerAuxiliary(iSubband + nSubbands) * ctranspose(powerCoefMatrix{iSubband + nSubbands})));
        end
        coefMatrix = coefMatrix + (3 / 2 * beta4) * powerRatio ^ 2 * (powerAuxiliary(nSubbands) * infoCoefMatrix{nSubbands} + infoAuxiliary(nSubbands) * powerCoefMatrix{nSubbands});

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression current;
            expression snr(nSubbands, 1);
            % \tilde(z)
            current = real(trace(coefMatrix * irsMatrix)) ...
                - (3 / 8 * beta4) * powerRatio ^ 2 * real(2 * (infoAuxiliary' * infoAuxiliary) + (powerAuxiliary' * powerAuxiliary)) ...
                - (3 / 2 * beta4) * powerRatio ^ 2 * real(infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));
            % R
            for iSubband = 1 : nSubbands
                snr(iSubband) = infoRatio * square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower;
            end
            maximize current;
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
                geo_mean(1 + snr) >= 2 ^ (rateConstraint / nSubbands);
        cvx_end

        % * Test convergence
        isConverged = (current - current_) / current <= tolerance || current == 0 || isnan(current);
        current_ = current;
    end

    % * Recover rank-1 solution by randomization method
    [u, sigma] = eig(irsMatrix);
    current = 0;
    for iCandidate = 1 : nCandidates
        irs_ = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
        irsMatrix = irs_ * irs_';
        % t_{I/P,n}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
        end
        % z
        current_ = real(1 / 2 * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + 3 / 8 * beta4 * powerRatio ^ 2 * (2 * (infoAuxiliary' * infoAuxiliary) + (powerAuxiliary' * powerAuxiliary)) ...
            * 3 / 2 * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands));
        % R
        rate_ = 0;
        for iSubband = 1 : nSubbands
            rate_ = rate_ + log2(1 + infoRatio * square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower);
        end
        if current_ > current && rate_ > rateConstraint
            current = current_;
            rate = rate_;
            irs = irs_;
        end
    end
    irs = irs(1 : nReflectors) / irs(end);

end
