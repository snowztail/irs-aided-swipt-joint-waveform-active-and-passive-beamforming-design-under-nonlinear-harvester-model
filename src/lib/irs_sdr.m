function [irsMatrix] = irs_sdr(k2, k4, resistance, noisePower, currentConstraint, tolerance, concatVector, concatMatrix, infoWaveform, powerWaveform, infoRatio, powerRatio)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %   - compute the output DC current and user rate
    %
    % Input:
    %   - k2: diode k-parameters
    %   - k4: diode k-parameters
    %   - resistance (R_ant): antenna resistance
    %   - noisePower (\sigma_n^2): average noise power
    %   - currentConstraint (z_0): average output DC current constraint
    %   - tolerance (\epsilon): minimum rate gain ratio per iteration
    %   - concatMatrix (R_n) [(nReflectors + 1) * (nReflectors + 1)]: rate SDR matrix
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (previous solution)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (previous solution)
    %   - infoRatio (\bar{\rho}): information splitting ratio (previous solution)
    %   - powerRatio (\rho): power splitting ratio (previous solution)
    %
    % Output:
    %   - irsMatrix: high-rank IRS outer product matrix
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    % * Construct current SDR matrices
    nSubbands = size(infoWaveform, 1);
    nReflectors = size(concatMatrix{1}, 1) - 1;
    % \boldsymbol{M}_{I/P}
    for iSubband = 1 : nSubbands
        infoMatrix = conj(infoWaveform) * transpose(infoWaveform);
        powerMatrix = conj(powerWaveform) * transpose(powerWaveform);
    end
    % \boldsymbol{C}_{I/P,n}
    infoSubmatrix = cell(2 * nSubbands - 1, 1);
    powerSubmatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoSubmatrix{iSubband + nSubbands} = concatVector' * diag(diag(infoMatrix, iSubband), iSubband) * concatVector;
        powerSubmatrix{iSubband + nSubbands} = concatVector' * diag(diag(powerMatrix, iSubband), iSubband) * concatVector;
    end

    % * SCA
    rate_ = 0;
    current_ = 0;
    isConverged = false;
    irsMatrix = eye(nReflectors + 1);
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    while ~isConverged
        % t_{I/P,n}^{(k-1)}
        for iSubband = - nSubbands + 1 : nSubbands - 1
            infoAuxiliary(iSubband + nSubbands) = trace(infoSubmatrix{iSubband + nSubbands} * irsMatrix);
            powerAuxiliary(iSubband + nSubbands) = trace(powerSubmatrix{iSubband + nSubbands} * irsMatrix);
        end
        % \boldsymbol{A}^{(k)}
        coefMatrix = (k2 * powerRatio * resistance / 2) * (infoSubmatrix{nSubbands} + powerSubmatrix{nSubbands});
        for iSubband = - nSubbands + 1 : nSubbands - 1
            coefMatrix = coefMatrix + (3 * k4 * powerRatio ^ 2 * resistance ^ 2 / 8) ...
                * (2 * (conj(infoAuxiliary(iSubband + nSubbands)) * infoSubmatrix{iSubband + nSubbands} + infoAuxiliary(iSubband + nSubbands) * ctranspose(infoSubmatrix{iSubband + nSubbands})) ...
                + (conj(powerAuxiliary(iSubband + nSubbands)) * powerSubmatrix{iSubband + nSubbands} + powerAuxiliary(iSubband + nSubbands) * ctranspose(powerSubmatrix{iSubband + nSubbands})));
        end
        coefMatrix = coefMatrix + (3 * k4 * powerRatio ^ 2 * resistance ^ 2 / 2) * (powerAuxiliary(nSubbands) * infoSubmatrix{nSubbands} + infoAuxiliary(nSubbands) * powerSubmatrix{nSubbands});
        % ensure strictly hermitian
        coefMatrix = (coefMatrix + coefMatrix') / 2;

        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression rate;
            expression current;
            % R
            for iSubband = 1 : nSubbands
                rate = rate + log(1 + infoRatio * square_abs(infoWaveform(iSubband)) * trace(concatMatrix{iSubband} * irsMatrix) / noisePower) / log(2);
            end
            % \tilde(z)
            current = trace(coefMatrix * irsMatrix);
            for iSubband = - nSubbands + 1 : nSubbands - 1
                current = current - (3 * k4 * powerRatio ^ 2 * resistance ^ 2 / 8) * (2 * infoAuxiliary(iSubband + nSubbands) * conj(infoAuxiliary(iSubband + nSubbands)) + powerAuxiliary(iSubband + nSubbands) * conj(powerAuxiliary(iSubband + nSubbands)));
            end
            current = current - (3 * k4 * powerRatio ^ 2 * resistance ^ 2 / 2) * real(infoAuxiliary(iSubband + nSubbands) * powerAuxiliary(iSubband + nSubbands));
            maximize rate;
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
                current >= currentConstraint;
        cvx_end

        % * Test convergence
        isConverged = (rate - rate_) / rate <= tolerance && ((current - current_) / current <= tolerance || current == 0);
        rate_ = rate
        current_ = current
    end

end
