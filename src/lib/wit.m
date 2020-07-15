function [capacity, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize user rate
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): transmit power constraint
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - capacity (R): channel capacity
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %   - construct information waveform by water-filling algorithm and matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % * Get data
    [nSubbands, nTxs, nReflectors] = size(incidentChannel);

    % * Initialize IRS and composite channel
    irs = ones(nReflectors, 1);
%     irs = exp(1i * 2 * pi * rand(nReflectors, 1));
    [compositeChannel, ~, concatSubchannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Construct waveform (water-filling + MRT) and splitting ratio
    [~, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
    infoWaveform = zeros(nTxs, nSubbands);
    for iSubband = 1 : nSubbands
        infoWaveform(:, iSubband) = sqrt(subbandPower(iSubband)) * compositeChannel(iSubband, :)' / norm(compositeChannel(iSubband, :));
    end
    powerWaveform = zeros(nTxs, nSubbands) + eps;
    infoRatio = 1 - eps;
    powerRatio = 1 - infoRatio;

    % * Obtain coefficients
    % \boldsymbol{M}_n, \boldsymbol{C}_n
    stackSubmatrix = cell(nSubbands, 1);
    rateMatrix = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        stackSubmatrix{iSubband} = [concatSubchannel{iSubband}; directChannel(iSubband, :)];
        rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
        rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
    end

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
            for iSubband = 1 : nSubbands
                snr(iSubband) = trace(rateMatrix{iSubband} * irsMatrix) / noisePower;
            end
            maximize geo_mean(1 + snr)
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
        cvx_end

        % * Recover rank-1 solution by randomization method
        [u, sigma] = eig(irsMatrix);
        rate = 0;
        for iCandidate = 1 : nCandidates
            irsCandidate = exp(1i * angle(u * sigma ^ (1 / 2) * (randn(nReflectors + 1, 1) + 1i * randn(nReflectors + 1, 1))));
            irsMatrix = irsCandidate * irsCandidate';

            % * Compute rate
            snr = zeros(nSubbands, 1);
            for iSubband = 1 : nSubbands
                snr(iSubband) = real(trace(rateMatrix{iSubband} * irsMatrix)) / noisePower;
            end
            rateCandidate = sum(log2(1 + snr));
            if rateCandidate > rate
                rate = rateCandidate;
                irs = irsCandidate;
            end
        end
        irs = irs(1 : nReflectors) / irs(end);

        % * Update composite channel and optimal waveform
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [capacity, subbandPower] = channel_capacity(compositeChannel, txPower, noisePower);
        for iSubband = 1 : nSubbands
            infoWaveform(:, iSubband) = sqrt(subbandPower(iSubband)) * compositeChannel(iSubband, :)' / norm(compositeChannel(iSubband, :));
        end

        % * Update coefficients
        for iSubband = 1 : nSubbands
            rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
            rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
        end

        % * Test convergence
        isConverged = abs(capacity - capacity_) / capacity <= tolerance;
        capacity_ = capacity;
    end

end
