function [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(channel, txPower)
    % Function:
    %   - initialize waveform and splitting ratio for GP algorithm
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): transmit power constraint
    %
    % Output:
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - the iterative GP algorithm is sensitive to initialization (only converges to stationary point)
    %   - use matched filter to initialize both waveform
    %   - set info splitting ratio to 1 as cvx crashes when initial rate is far from constraint
    %   - GP requires nonzero entries thus use eps to replace zero
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Jul 20


    % * Get data
    nSubbands = size(channel, 1);

    % * Initialize algorithm
    infoRatio = 1 - eps;
    powerRatio = eps;
    infoWaveform = sqrt(txPower / nSubbands) * channel' ./ vecnorm(channel, 2, 2)';
    powerWaveform = sqrt(txPower / nSubbands) * channel' ./ vecnorm(channel, 2, 2)';

end
