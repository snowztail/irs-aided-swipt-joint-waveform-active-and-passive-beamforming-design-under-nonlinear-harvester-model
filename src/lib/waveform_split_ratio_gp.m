function [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_split_ratio_gp(beta2, beta4, channel, infoWaveform, powerWaveform, infoRatio, powerRatio, txPower, noisePower, rateConstraint, tolerance)
    % Function:
    %   - jointly optimize waveform and splitting ratio to maximize the R-E region
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform (in the previous iteration)
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform (in the previous iteration)
    %   - infoRatio (\bar{\rho}): information splitting ratio (in the previous iteration)
    %   - powerRatio (\rho): power splitting ratio (in the previous iteration)
    %   - txPower (P): transmit power constraint
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - the optimal spatial single-user beamformer is MRT
    %   - effective channel is given by Euclidean norm
    %   - only power allocation (amplitude optimization) problem in frequency domain (DoF = N)
    %   - AM-GM is only suitable for real variables, thus not for multiuser case (need to optimize spatial beamformer)
    %   - tightness of the AM-GM equality depends on exponentials
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % * Get data
    [nSubbands] = size(channel, 1);

    % * Initialize algorithm
    % \boldsymbol{a}, \boldsymbol{s}_{I/P}
    channelAmplitude = vecnorm(channel, 2, 2);
    infoAmplitude = zeros(nSubbands, 1);
    powerAmplitude = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        infoAmplitude(iSubband) = abs(channel(iSubband, :) / norm(channel(iSubband, :)) * infoWaveform(:, iSubband));
        powerAmplitude(iSubband) = abs(channel(iSubband, :) / norm(channel(iSubband, :)) * powerWaveform(:, iSubband));
    end
    % \gamma_{I/P}
    [~, ~, currentExponent] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
    [~, ~, rateExponent] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

    % * Iterative GP
    isConverged = false;
    current_ = 0;
    while ~isConverged
        cvx_begin gp quiet
            cvx_solver mosek
            variable auxiliary
            variable infoAmplitude(nSubbands, 1) nonnegative
            variable powerAmplitude(nSubbands, 1) nonnegative
            variable infoRatio nonnegative
            variable powerRatio nonnegative

            [~, currentMonomial, ~] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
            [~, rateMonomial, ~] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

            minimize (1 / auxiliary)
            subject to
                (1 / 2) * (infoAmplitude' * infoAmplitude + powerAmplitude' * powerAmplitude) <= txPower;
                auxiliary * prod((currentMonomial ./ currentExponent) .^ (-currentExponent)) <= 1;
                2 ^ rateConstraint * prod(prod((rateMonomial ./ rateExponent) .^ (-rateExponent))) <= 1;
                powerRatio + infoRatio <= 1;
        cvx_end

        [current, ~, currentExponent] = current_gp(beta2, beta4, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
        [rate, ~, rateExponent] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower);

        % * Test convergence
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end

    % * Reconstruct waveform by power allocation + beamforming
    infoWaveform = transpose(infoAmplitude) .* channel' ./ vecnorm(channel, 2, 2)';
    powerWaveform = transpose(powerAmplitude) .* channel' ./ vecnorm(channel, 2, 2)';

end
