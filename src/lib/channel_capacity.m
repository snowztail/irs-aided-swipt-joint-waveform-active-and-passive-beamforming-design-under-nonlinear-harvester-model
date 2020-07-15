function [capacity, subbandPower] = channel_capacity(channel, txPower, noisePower)
    % Function:
    %   - calculate the maximum achievable rate based on water-filling power allocation
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: the wireless channel
    %   - txPower (P): average transmit power
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity: the maximum achievable rate
    %   - subbandPower: optimal power allocation
    %
    % Comment:
    %   - for MISO OFDM channels
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Aug 19



    % * MRT within each subband
    subbandStrength = vecnorm(channel, 2, 2) .^ 2;

    % * Power allocation
    subbandPower = water_filling(subbandStrength, 2 * txPower, noisePower);

    % * Compute capacity
    capacity = sum(log2(1 + subbandPower .* subbandStrength / noisePower));

end
