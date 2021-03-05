function [sample, solution] = re_sample_swipt_gp(alpha, beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance)
    % Function:
    %   - sample R-E region by computing the output DC current and rate
    %
    % Input:
    %   - alpha: scale ratio of SMF
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - nSamples (S): number of samples in R-E region
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - AO algorithm only converge to stationary points
    %   - proceed from high rate points to high current points
    %   - results are sensitive to initialization
    %   - under default initialization, some samples may be strictly worse than previous ones (especially for a large number of transmit antennas or reflectors)
    %   - if the issue above happens, we discard the result based on default initialization and reinitialize this point by previous solution
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    sample = zeros(2, nSamples);
    solution = cell(nSamples, 1);

    % * Initialize algorithm by WIT point and set rate constraints
    [sample(:, 1), solution{1}] = re_sample_wit_wf(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
	irs = solution{1}.irs;
    compositeChannel = solution{1}.compositeChannel;
    rateConstraint = linspace(sample(1, 1), 0, nSamples);

    % * Non-WIT points
    for iSample = 2 : nSamples
        isDominated = false;
        while true
            if ~isDominated
                % * Default initialization
                [infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize_waveform(alpha, beta2, beta4, compositeChannel, txPower, noisePower);
            else
                % * Initialize with previous solution
                struct2variables(solution{iSample - 1});
            end
            [rate, current, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = waveform_gp(beta2, beta4, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), tolerance);
            [infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);

            % * Alternating optimization
            isConverged = false;
			current_ = 0;
			eigRatio = [];
            while ~isConverged
                [irs, eigRatio(end + 1)] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint(iSample), nCandidates, tolerance);
                [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
                [rate, current, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = waveform_gp(beta2, beta4, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), tolerance);
                [infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);
                isConverged = abs(current - current_) <= tolerance;
                current_ = current;
            end

            % * Check whether strictly dominated
            isDominated = current <= sample(2, iSample - 1);
            if ~isDominated
                break;
            end
        end
        sample(:, iSample) = [rate; current];
        solution{iSample} = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);
    end

end
