function [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(channel, txPower, noisePower)
    % Function:
    %   - initialize waveform and splitting ratio for GP algorithm
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - the iterative GP algorithm is sensitive to initialization (only converges to stationary point)
    %   - GP requires nonzero entries thus use eps to replace zero
    %   - cvx crashes when initial rate is far from constraint
    %       - use matched filter to initialize power waveform with power P
    %       - use water-filling + MRT to initialize information waveform with power P
    %       - initialize information splitting ratio to 1
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Jul 20


    % * Get data
    nSubbands = size(channel, 1);

    % * Initialize algorithm
    infoRatio = 1;
    powerRatio = 1;
    [~, infoWaveform] = channel_capacity(channel, txPower, noisePower);
    % infoWaveform = sqrt(txPower / nSubbands) * channel' ./ vecnorm(channel, 2, 2)';
    powerWaveform = sqrt(txPower / nSubbands) * channel' ./ vecnorm(channel, 2, 2)';

end
