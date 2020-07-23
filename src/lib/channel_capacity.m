function [capacity, infoWaveform] = channel_capacity(channel, txPower, noisePower)
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
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %
    % Comment:
    %   - for MISO OFDM channels
    %   - multiply 2 in power allocation to convert effective value to magnitude
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Aug 19



    % * MRT within each subband
    subbandStrength = vecnorm(channel, 2, 2) .^ 2;

    % * Water-filling power allocation
    subbandPower = water_filling(subbandStrength, 2 * txPower, noisePower);

    % * Construct waveform
    infoWaveform = sqrt(subbandPower') .* channel' ./ vecnorm(channel, 2, 2)';

    % * Compute capacity
    capacity = sum(log2(1 + subbandPower .* subbandStrength / noisePower));

end
