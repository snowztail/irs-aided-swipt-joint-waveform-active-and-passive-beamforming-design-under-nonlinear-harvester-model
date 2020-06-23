function [infoRatio, powerRatio] = split_ratio(channel, noisePower, rateConstraint, infoWaveform)
    % Function:
    %   - optimize the information and power splitting ratio to maximize the R-E region
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): user rate constraint
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %
    % Output:
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - tight the rate constraint
    %
    % Author & Date: Yang (i@snowztail.com) - 20 Jun 20



    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);

    % * Solve splitting ratio by CVX
    cvx_begin quiet
        cvx_solver mosek
        cvx_precision high
        variable powerRatio nonnegative
        expression snr(nSubbands, 1);
        % \gamma
        for iSubband = 1 : nSubbands
            snr(iSubband) = (1 - powerRatio) * abs(channel(iSubband) * infoWaveform(iSubband)) ^ 2 / noisePower;
        end
        maximize powerRatio
        subject to
            geo_mean(1 + snr) >= 2 ^ (rateConstraint / nSubbands);
            powerRatio <= 1;
    cvx_end
    infoRatio = 1 - powerRatio;

end
