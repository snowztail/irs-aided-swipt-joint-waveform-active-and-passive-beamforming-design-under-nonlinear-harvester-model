function [sample, solution] = re_sample_fixed_irs(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, nSamples, tolerance)
    % Function:
    %   - sample R-E region with fixed WIT-optimized IRS by computing the output DC current and rate
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - nSamples (S): number of samples in R-E region
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: waveform and splitting ratio
    %
    % Comment:
    %   - fix IRS, only optimize waveform
    %   - suboptimal algorithm only converge to stationary points
    %   - proceed from high rate points to high current points
    %   - results are sensitive to initialization
    %   - under default initialization, some samples may be strictly worse than previous ones (especially for a large number of transmit antennas or reflectors)
    %   - if the issue above happens, we discard the result based on default initialization and reinitialize this point by previous solution
    %
    % Author & Date: Yang (i@snowztail.com) - 7 Aug 20


    sample = zeros(2, nSamples);
    solution = cell(nSamples, 1);

    % * Initialize algorithm and set rate constraints
    [capacity, irs, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = wit(directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance);
    [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    rateConstraint = linspace(capacity, 0, nSamples);

    % * WIT point
    sample(:, 1) = [capacity; 0];
    solution{1} = variables2struct(infoAmplitude, powerAmplitude, infoRatio, powerRatio);

    % * Non-WIT points
    for iSample = 2 : nSamples
        isDominated = false;
        while true
            if ~isDominated
                % * Default initialization
                [infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize_waveform(compositeChannel, txPower, noisePower);
            else
                % * Initialize with previous solution
                struct2variables(solution{iSample - 1});
            end
            % * Optimize waveform with WIT-optimized IRS
            [infoAmplitude, powerAmplitude, infoRatio, powerRatio, rate, current] = waveform_gp(beta2, beta4, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, txPower, noisePower, rateConstraint(iSample), tolerance);

            % * Check whether strictly dominated
            isDominated = current <= sample(2, iSample - 1);
            if ~isDominated
                break;
            end
        end
        sample(:, iSample) = [rate; current];
        solution{iSample} = variables2struct(infoAmplitude, powerAmplitude, infoRatio, powerRatio);
    end

end
