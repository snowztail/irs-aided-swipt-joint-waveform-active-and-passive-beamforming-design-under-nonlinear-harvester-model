function [irs, eigRatio, scaIter] = irs_sdr(beta2, beta4, directChannel, cascadedChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates, tolerance)
    % Function:
    %   - optimize the IRS reflection coefficients to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
    %   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %   - irs (\phi) [nReflectors * 1]: IRS reflection coefficients (in the previous iteration)
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
	%   - irs (\phi) [nReflectors * 1]: IRS reflection coefficients
	%	- eigRatio (r): the maximum eigenvalue of the relaxed solution over the sum eigenvalue of the relaxed solution
	%	- scaIter: R-E sample at each iteration
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20


    % * Get data
    [nReflectors, nTxs, nSubbands] = size(cascadedChannel);

    % * Construct coefficient matrices
    % \boldsymbol{W}_{I/P}, \boldsymbol{W}_{I/P,n}
    infoMatrix = vec(infoWaveform) * vec(infoWaveform)';
    powerMatrix = vec(powerWaveform) * vec(powerWaveform)';
    infoBlkDiag = block_diagonal(infoMatrix, nTxs, nSubbands);
    powerBlkDiag = block_diagonal(powerMatrix, nTxs, nSubbands);

    % \boldsymbol{M}_n, \boldsymbol{C}_n
    stackSubmatrix = cell(nSubbands, 1);
    rateMatrix = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
		stackSubmatrix{iSubband} = [cascadedChannel(:, :, iSubband); directChannel(iSubband, :)];
        rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
        rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
    end

    % \boldsymbol{M}, \boldsymbol{C}_{I/P,n}
    stackMatrix = [reshape(cascadedChannel, [nReflectors, nTxs * nSubbands]); vec(directChannel')'];
    infoCurrentMatrix = cell(2 * nSubbands - 1, 1);
    powerCurrentMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoCurrentMatrix{iSubband + nSubbands} = stackMatrix * infoBlkDiag{iSubband + nSubbands} * stackMatrix';
        powerCurrentMatrix{iSubband + nSubbands} = stackMatrix * powerBlkDiag{iSubband + nSubbands} * stackMatrix';
    end
    infoCurrentMatrix{nSubbands} = hermitianize(infoCurrentMatrix{nSubbands});
    powerCurrentMatrix{nSubbands} = hermitianize(powerCurrentMatrix{nSubbands});

    % * Initialize IRS
    % \bar{\boldsymbol{\phi}}, bar{\boldsymbol{\Phi}}
    irs = [irs; 1];
    irsMatrix = irs * irs';

    % * Initialize auxiliaries
    % t_{I/P,n}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(infoCurrentMatrix{iSubband + nSubbands} * irsMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(powerCurrentMatrix{iSubband + nSubbands} * irsMatrix);
    end
    infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
    powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

	% * Initial R-E point
	current_ = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);
	snr = zeros(nSubbands, 1);
	for iSubband = 1 : nSubbands
		snr(iSubband) = infoRatio * real(trace(rateMatrix{iSubband} * irsMatrix)) / noisePower;
	end
	rate_ = sum(log2(1 + snr));
	scaIter = [rate_; current_];

    % * SCA
    isConverged = false;
    while ~isConverged
        % * Update auxiliaries
        infoAuxiliary_ = infoAuxiliary;
        powerAuxiliary_ = powerAuxiliary;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            cvx_precision high
            cvx_expert true
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression infoAuxiliary(2 * nSubbands - 1, 1);
            expression powerAuxiliary(2 * nSubbands - 1, 1);
            expression snr(nSubbands, 1);
            % t_{I/P,n}
            for iSubband = - nSubbands + 1 : nSubbands - 1
                infoAuxiliary(iSubband + nSubbands) = trace(infoCurrentMatrix{iSubband + nSubbands} * irsMatrix);
                powerAuxiliary(iSubband + nSubbands) = trace(powerCurrentMatrix{iSubband + nSubbands} * irsMatrix);
            end
            infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
            powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));
            % \tilde{z}
            currentSca = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
                + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * (2 * infoAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) ^ 2) + 2 * real(powerAuxiliary_' * powerAuxiliary) - powerAuxiliary_' * powerAuxiliary_) ...
                + (3 / 2) * beta4 * powerRatio ^ 2 * (infoAuxiliary(nSubbands) * powerAuxiliary_(nSubbands) + powerAuxiliary(nSubbands) * infoAuxiliary_(nSubbands) - infoAuxiliary_(nSubbands) * powerAuxiliary_(nSubbands));
            % \gamma
            for iSubband = 1 : nSubbands
                snr(iSubband) = infoRatio * trace(rateMatrix{iSubband} * irsMatrix) / noisePower;
            end
            % R
            rate = sum(log(1 + snr) / log(2));
            maximize currentSca;
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
                rate >= rateConstraint;
        cvx_end

        % * Update output current
        current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Test convergence
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
		scaIter(:, end + 1) = [rate; current];
    end
    irsMatrix = full(irsMatrix);

	% * Examine the unit-rank property of the result by eigenvalue ratio
	[u, sigma] = eig(irsMatrix);
	eigRatio = max(sigma) / sum(sigma);

	% * Recover rank-1 solution by randomization method
    current = 0;
    for iCandidate = 1 : nCandidates
        irsCandidate = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
        irsMatrix = irsCandidate * irsCandidate';

        % * Update auxiliaries
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(infoCurrentMatrix{iSubband + nSubbands} * irsMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerCurrentMatrix{iSubband + nSubbands} * irsMatrix);
        end
        infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
        powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

        % * Compute rate and current
        for iSubband = 1 : nSubbands
            snr(iSubband) = infoRatio * real(trace(rateMatrix{iSubband} * irsMatrix)) / noisePower;
        end
        rateCandidate = sum(log2(1 + snr));

        currentCandidate = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
            + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
            + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

        % * Choose best candidate
        if currentCandidate >= current && rateCandidate >= rateConstraint
            current = currentCandidate;
            irs = irsCandidate;
        end
    end
    irs = reshape(irs(1 : nReflectors) / irs(end), [nReflectors, 1]);

end
