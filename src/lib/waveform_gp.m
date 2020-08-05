function [infoAmplitude, powerAmplitude, infoRatio, powerRatio, rate, current] = waveform_gp(beta2, beta4, channel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, txPower, noisePower, rateConstraint, tolerance)
    % Function:
    %   - jointly optimize waveform and splitting ratio to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain (in the previous iteration)
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain (in the previous iteration)
    %   - infoRatio (\bar{\rho}): information splitting ratio (in the previous iteration)
    %   - powerRatio (\rho): power splitting ratio (in the previous iteration)
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - obtain waveform amplitude in frequency domain
    %   - effective channel is given by Euclidean norm
    %   - only power allocation (amplitude optimization) problem in frequency domain (DoF = N)
    %   - AM-GM is only suitable for real variables, thus not for multiuser case (need to optimize spatial beamformer)
    %   - tightness of the AM-GM equality depends on exponentials
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    % * Get data
    nSubbands = size(channel, 1);

    % * Initialize algorithm
    % \boldsymbol{a}
    channelAmplitude = vecnorm(channel, 2, 2);
    % \gamma_{I/P}
    [~, ~, currentExponent] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
    [~, ~, rateExponent] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

    % * Iterative GP
    isConverged = false;
    current_ = 0;
    while ~isConverged
        cvx_begin gp quiet
            cvx_solver mosek
            % cvx_precision high
            variable auxiliary
            variable infoAmplitude(1, nSubbands) nonnegative
            variable powerAmplitude(1, nSubbands) nonnegative
            variable infoRatio nonnegative
            variable powerRatio nonnegative

            [~, currentMonomial, ~] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
            [~, rateMonomial, ~] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

            minimize (1 / auxiliary)
            subject to
                (1 / 2) * (norm(infoAmplitude) ^ 2 + norm(powerAmplitude) ^ 2) <= txPower;
                auxiliary * prod((currentMonomial ./ currentExponent) .^ (-currentExponent)) <= 1;
                2 ^ rateConstraint * prod(prod((rateMonomial ./ rateExponent) .^ (-rateExponent))) <= 1;
                powerRatio + infoRatio <= 1;
        cvx_end

        [current, ~, currentExponent] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
        [rate, ~, rateExponent] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

        % * Test convergence
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
    end

end
