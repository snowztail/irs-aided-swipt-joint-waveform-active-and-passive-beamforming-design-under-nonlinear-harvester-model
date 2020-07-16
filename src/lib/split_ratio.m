function [infoRatio, powerRatio] = split_ratio(channel, infoWaveform, noisePower, rateConstraint)
    % Function:
    %   - optimize the information and power splitting ratio to maximize the R-E region
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %
    % Output:
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - tight the rate constraint
    %
    % Author & Date: Yang (i@snowztail.com) - 20 Jun 20



    % * Get data
    nSubbands = size(channel, 1);

    % * Solve splitting ratio by CVX
    cvx_begin quiet
        cvx_solver mosek
        cvx_precision high
        variable powerRatio nonnegative
        expression snr(nSubbands, 1);
        infoRatio = 1 - powerRatio;
        % \gamma
        for iSubband = 1 : nSubbands
            snr(iSubband) = infoRatio * abs(channel(iSubband, :) * infoWaveform(:, iSubband)) ^ 2 / noisePower;
        end
        maximize powerRatio
        subject to
            geo_mean(1 + snr) >= 2 ^ (rateConstraint / nSubbands);
            powerRatio <= 1;
    cvx_end

end
