function [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, txPower, rateConstraint, tolerance, infoRatio, powerRatio, noisePower, channel, infoWaveform, powerWaveform)
    % Function:
    %   - jointly optimize waveform and splitting ratio to maximize the R-E region
    %   - based on GP
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - txPower (P): transmit power constraint
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - noisePower (\sigma_n^2): average noise power
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (in the previous iteration)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (in the previous iteration)
    %   - infoRatio (\bar{\rho}): information splitting ratio (in the previous iteration)
    %   - powerRatio (\rho): power splitting ratio (in the previous iteration)
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - rate (R): user sum rate over all subbands
    %   - current (z): harvester output DC current
    %
    % Comment:
    %   - the tightness of the AM-GM equality depends on the exponentials
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);
    % \psi_{I/P,n}
    phase = -angle(channel);
    % \boldsymbol{A}
    channel = abs(channel);
    % \boldsymbol{s}_{I/P}
    infoWaveform = abs(infoWaveform);
    powerWaveform = abs(powerWaveform);

    % * Iterative GP
    isConverged = false;
    current_ = 0;
    [~, ~, currentExponent] = current_gp(beta2, beta4, channel, infoWaveform, powerWaveform, powerRatio);
    [~, ~, rateExponent] = rate_gp(noisePower, channel, infoWaveform, infoRatio);
    while ~isConverged
        cvx_begin gp quiet
            cvx_solver mosek
            variable auxiliary
            variable infoWaveform(nSubbands, 1) nonnegative
            variable powerWaveform(nSubbands, 1) nonnegative
            variable infoRatio nonnegative
            variable powerRatio nonnegative

            [~, currentMonomial, ~] = current_gp(beta2, beta4, channel, infoWaveform, powerWaveform, powerRatio);
            [~, rateMonomial, ~] = rate_gp(noisePower, channel, infoWaveform, infoRatio);

            minimize (1 / auxiliary)
            subject to
                (1 / 2) * (infoWaveform' * infoWaveform + powerWaveform' * powerWaveform) <= txPower;
                auxiliary * prod((currentMonomial ./ currentExponent) .^ (-currentExponent)) <= 1;
                2 ^ rateConstraint * prod(prod((rateMonomial ./ rateExponent) .^ (-rateExponent))) <= 1;
                powerRatio + infoRatio <= 1;
        cvx_end

        [current, ~, currentExponent] = current_gp(beta2, beta4, channel, infoWaveform, powerWaveform, powerRatio);
        [rate, ~, rateExponent] = rate_gp(noisePower, channel, infoWaveform, infoRatio);

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance || current == 0;
        current_ = current;
    end
    infoWaveform = infoWaveform .* exp(1i * phase);
    powerWaveform = powerWaveform .* exp(1i * phase);

end
