function [capacity, subbandPower] = channel_capacity(channel, txPower, noisePower)
    % Function:
    %   - calculate the maximum achievable rate based on water-filling power allocation
    %
    % Input:
    %   - channel (h) [nSubbands]: the wireless channel
    %   - txPower (P): average transmit power
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity: the maximum achievable rate
    %   - subbandPower: optimal power allocation
    %
    % Comment:
    %   - for SISO OFDM channels
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Aug 19



    subbandStrength = abs(channel) .^ 2;
    subbandPower = water_filling(subbandStrength, 2 * txPower, noisePower);
    capacity = sum(log2(1 + subbandPower .* subbandStrength / noisePower));

end
