function [sample, solution] = re_sample_wit_wf(directChannel, cascadedChannel, txPower, noisePower, nCandidates, tolerance)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize user rate
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
	%   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, composite channel, waveform, splitting ratio, eigenvalue ratio, R-E sample at each iteration
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %   - construct information waveform by water-filling algorithm and matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    % * Get data
    [nReflectors, ~, nSubbands] = size(cascadedChannel);

    % * Initialize IRS and composite channel
    irs = exp(1i * 2 * pi * rand(nReflectors, 1));
    [compositeChannel] = composite_channel(directChannel, cascadedChannel, irs);

    % * Construct waveform (water-filling + MRT) and initialize splitting ratio
    [capacity_, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = water_filling(compositeChannel, txPower, noisePower);
    [infoWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);

	% * Initial R-E point
	bcdIter = [capacity_; 0];

    % * Obtain coefficients
    % \boldsymbol{M}_n, \boldsymbol{C}_n
    stackSubmatrix = cell(nSubbands, 1);
    rateMatrix = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        stackSubmatrix{iSubband} = [cascadedChannel(:, :, iSubband); directChannel(iSubband, :)];
        rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
        rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
    end

    % * AO
    isConverged = false;
	eigRatio = [];
    while ~isConverged
        % * Solve high-rank outer product matrix by CVX
        cvx_begin quiet
            cvx_solver mosek
            cvx_precision high
            cvx_expert true
            variable irsMatrix(nReflectors + 1, nReflectors + 1) hermitian semidefinite;
            expression snr(nSubbands, 1);
            for iSubband = 1 : nSubbands
                snr(iSubband) = trace(rateMatrix{iSubband} * irsMatrix) / noisePower;
            end
            rate = sum(log(1 + snr) / log(2));
            maximize rate
            subject to
                diag(irsMatrix) == ones(nReflectors + 1, 1);
        cvx_end
        irsMatrix = full(irsMatrix);

		% * Examine the unit-rank property of the result by eigenvalue ratio
		[u, sigma] = eig(irsMatrix);
		eigRatio(end + 1) = max(sigma) / sum(sigma);

		% * Recover rank-1 solution by randomization method
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

            % * Choose best candidate
            if rateCandidate > rate
                rate = rateCandidate;
                irs = irsCandidate;
            end
        end
        irs = reshape(irs(1 : nReflectors) / irs(end), [nReflectors, 1]);

        % * Update composite channel and optimal waveform
        [compositeChannel] = composite_channel(directChannel, cascadedChannel, irs);
        [capacity, infoAmplitude] = water_filling(compositeChannel, txPower, noisePower);
        [infoWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);

        % * Update coefficients
        for iSubband = 1 : nSubbands
            rateMatrix{iSubband} = stackSubmatrix{iSubband} * infoWaveform(:, iSubband) * infoWaveform(:, iSubband)' * stackSubmatrix{iSubband}';
            rateMatrix{iSubband} = hermitianize(rateMatrix{iSubband});
        end

        % * Test convergence
        isConverged = abs(capacity - capacity_) <= tolerance;
        capacity_ = capacity;
		bcdIter(:, end + 1) = [capacity; 0];
    end

	sample = [capacity; eps];
	solution = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio, bcdIter);

end
