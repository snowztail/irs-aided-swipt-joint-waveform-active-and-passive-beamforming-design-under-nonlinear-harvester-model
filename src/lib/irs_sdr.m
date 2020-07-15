function [irs] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates, tolerance)
    % Function:
    %   - optimize the IRS reflection coefficients to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients (in the previous iteration)
    %
    % Output:
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20


    % * Get data
    [nSubbands, nTxs, nReflectors] = size(incidentChannel);
    [compositeChannel, concatChannel, concatSubchannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Construct coefficient matrices
    % \boldsymbol{W}_{I/P}
    infoMatrix = vec(infoWaveform) * vec(infoWaveform)';
    powerMatrix = vec(powerWaveform) * vec(powerWaveform)';

    iM_ = cell(nSubbands, 1);
    pM_ = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        iM_{iSubband} = zeros(nTxs * nSubbands);
        pM_{iSubband} = zeros(nTxs * nSubbands);
        for jSubband = 1 : nSubbands + 1 - iSubband
            iM_{iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = infoMatrix((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            pM_{iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = powerMatrix((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
        end
    end

    iM = cell(2 * nSubbands - 1, 1);
    pM = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        if iSubband >= 0
            iM{iSubband + nSubbands} = iM_{iSubband + 1};
            pM{iSubband + nSubbands} = pM_{iSubband + 1};
        else
            iM{iSubband + nSubbands} = iM_{-iSubband + 1}';
            pM{iSubband + nSubbands} = pM_{-iSubband + 1}';
        end
    end

    % \boldsymbol{M}_n, \boldsymbol{C}_n
    stackSubmatrix = cell(nSubbands, 1);
    rateMatrix = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        stackSubmatrix{iSubband} = [concatSubchannel{iSubband}; directChannel(iSubband, :)];
        rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
        rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
    end

    % \boldsymbol{M}, \boldsymbol{C}_{I/P,n}
    stackMatrix = [concatChannel; vec(directChannel')'];
    infoCurrentMatrix = cell(2 * nSubbands - 1, 1);
    powerCurrentMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoCurrentMatrix{iSubband} = stackMatrix * conj(iM{iSubband + nSubbands}) * stackMatrix';
        powerCurrentMatrix{iSubband} = stackMatrix * conj(pM{iSubband + nSubbands}) * stackMatrix';
    end




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
    infoCoefMatrix{nSubbands} = hermitianize(infoCoefMatrix{nSubbands});
    powerCoefMatrix{nSubbands} = hermitianize(powerCoefMatrix{nSubbands});
    % \bar{\boldsymbol{\phi}}^{(0)}
    irs = [irs; 1];
    % \bar{\boldsymbol{\Phi}}^{(0)}
    irsMatrix = irs * irs';
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
        % * Update auxiliary variables
        % t_{I/P,n}^{(i-1)}
        infoAuxiliary_ = infoAuxiliary;
        powerAuxiliary_ = powerAuxiliary;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            cvx_precision high
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            expression snr(nSubbands, 1);
            % t_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(infoCoefMatrix{iSubband + nSubbands} * irsMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(powerCoefMatrix{iSubband + nSubbands} * irsMatrix);
            end
            infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
            powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));
            % \tilde{z}
            currentLb = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + 2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_) ...
                + (3 / 2) * beta4 * powerRatio ^ 2 * ((1 / 2) * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) - (1 / 4) * (infoAuxiliary_(nSubbands) + powerAuxiliary_(nSubbands)) ^ 2 - (1 / 4) * (infoAuxiliary(nSubbands) - powerAuxiliary(nSubbands)) ^ 2);
            % \gamma
            for iSubband = 1 : nSubbands
                snr(iSubband) = infoRatio * abs(infoWaveform(iSubband)) ^ 2 * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower;
            end
%             rate = sum_log(1 + snr) / log(2);
            maximize currentLb;
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
                geo_mean(1 + snr) >= 2 ^ (rateConstraint / nSubbands);
%                 rate >= rateConstraint;
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
    irsMatrix = full(irsMatrix);

    % * Recover rank-1 solution by randomization method
    [u, sigma] = eig(irsMatrix);
    current = 0;
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
        infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
        powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));
        % z
        current_ = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);
        % \bar{\boldsymbol{\phi}}^\star
        if current_ >= current
            current = current_;
            irs = irs_;
        end
    end
    % \boldsymbol{\phi}
    irs = irs(1 : nReflectors) / irs(end);

end
