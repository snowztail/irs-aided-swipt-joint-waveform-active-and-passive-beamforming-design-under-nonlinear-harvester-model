function [irs, rate] = irs_flat_wit(noisePower, concatMatrix, infoWaveform, nCandidates)
    % Function:
    %   - optimize the IRS reflection coefficients to maximize user rate
    %
    % Input:
    %   - noisePower (\sigma_n^2): average noise power
    %   - concatMatrix (R_n) [(nReflectors + 1) * (nReflectors + 1)]: rate SDR matrix
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (in the previous iteration)
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %
    % Output:
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - rate: maximum achievable user rate
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %
    % Author & Date: Yang (i@snowztail.com) - 5 Jun 20



    % * Construct current SDR matrices
    nSubbands = size(infoWaveform, 1);
    nReflectors = size(concatMatrix{1}, 1) - 1;
    irsMatrix_ = ones(nReflectors + 1);
    compositeChannel_ = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        compositeChannel_(iSubband) = trace(concatMatrix{iSubband} * irsMatrix_);
    end
    cvx_begin quiet
        cvx_solver mosek
        variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
        expression compositeChannel(nSubbands, 1)
        expression rate;
        for iSubband = 1 : nSubbands
            compositeChannel(iSubband) = trace(concatMatrix{iSubband} * irsMatrix);
            rate = rate + log(1 + square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower) / log(2);
        end
        maximize rate
        subject to
            diag(irsMatrix) == ones(nReflectors + 1, 1);
    cvx_end
    irsMatrix = full(irsMatrix);

    % * Recover rank-1 solution by randomization method
    [u, sigma] = eig(irsMatrix);
    rate = 0;
    for iCandidate = 1 : nCandidates
        irs_ = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
        irsMatrix = irs_ * irs_';
        % R
        rate_ = 0;
        for iSubband = 1 : nSubbands
            rate_ = rate_ + log2(1 + square_abs(infoWaveform(iSubband)) * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower);
        end
        if rate_ > rate
            rate = rate_;
            irs = irs_;
        end
    end
    irs = irs(1 : nReflectors) / irs(end);

end
