function [capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_ff(irsGain, tolerance, directChannel, incidentChannel, reflectiveChannel, txPower, nCandidates, noisePower)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize user rate
    %
    % Input:
    %   - irsGain (\beta): gain on each reflecting element
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): transmit power constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity (R): channel capacity
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %   - construct information waveform by water-filling algorithm and matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % * Initialize IRS and waveform
    [nSubbands, nReflectors, ~] = size(reflectiveChannel);
    % \boldsymbol{\phi}^{(0)}
    irs = irsGain * ones(nReflectors, 1);
    % h^{(0)}, \boldsymbol{R}_n
    [compositeChannel, ~, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    % p_n^{(0)}
    [~, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    % \rho
    powerRatio = eps;
    % \bar{\rho}
    infoRatio = 1 - powerRatio;
    % \boldsymbol{W}_I^{(0)}
    infoWaveform = sqrt(subbandPower) .* exp(1i * angle(conj(compositeChannel)));
    % \boldsymbol{W}_P
    powerWaveform = zeros(size(compositeChannel)) + eps;

    % * AO
    isConverged = false;
    capacity_ = 0;
    while ~isConverged
        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            cvx_precision high
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression snr(nSubbands, 1);
            % \gamma
            for iSubband = 1 : nSubbands
                snr(iSubband) = abs(infoWaveform(iSubband)) ^ 2 * trace(concatMatrix{iSubband} * irsMatrix) / noisePower;
            end
            maximize geo_mean(1 + snr)
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
        cvx_end

        % * Recover rank-1 solution by randomization method
        [u, sigma] = eig(irsMatrix);
        rate = 0;
        for iCandidate = 1 : nCandidates
            irs_ = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
            irsMatrix = irs_ * irs_';
            % R
            rate_ = 0;
            for iSubband = 1 : nSubbands
                rate_ = rate_ + log2(1 + abs(infoWaveform(iSubband)) ^ 2 * real(trace(concatMatrix{iSubband} * irsMatrix)) / noisePower);
            end
            % \bar{\boldsymbol{\phi}}^\star
            if rate_ > rate
                rate = rate_;
                irs = irs_;
            end
        end
        irs = irs(1 : nReflectors) / irs(end);

        % * Construct information waveform
        % h^{(i)}
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        % p_n^{(i)}
        [capacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
        % \boldsymbol{W}_I^{(i)}
        infoWaveform = sqrt(subbandPower) .* exp(1i * angle(conj(compositeChannel)));

        % * Test convergence
        isConverged = abs(capacity - capacity_) / capacity <= tolerance;
        capacity_ = capacity;
    end

end
