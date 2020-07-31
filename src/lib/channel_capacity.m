function [capacity, infoWaveform] = channel_capacity(channel, txPower, noisePower)
    % Function:
    %   - calculate the maximum achievable rate based on water-filling power allocation
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity: the maximum achievable rate
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %
    % Comment:
    %   - for MISO OFDM channels
    %   - multiply 2 in power allocation to convert average power to maximum power
    %   - this is because wireless communication always transmit at the peak of waveform (maybe a different definition of transmit power budget)
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Aug 19


    % * Get data
    nSubbands = size(channel, 1);

    % * Convert average power to maximum power (different definition, transmit at peaks)
    txPower = 2 * txPower;
    snr = txPower / noisePower;

    % * Obtain subchannel strength (MRT within each subband)
    strength = vecnorm(channel, 2, 2) .^ 2;

    % * Iterative water-filling power allocation
    subbandPower = zeros(nSubbands, 1);
    [strength, index] = sort(strength, 'descend');
    baseLevel = 1 ./ (snr * strength);
    for iSubband = 1 : nSubbands
        waterLevel = 1 / (nSubbands - iSubband + 1) * (1 + sum(baseLevel(1 : (nSubbands - iSubband + 1))));
        candidate = waterLevel - baseLevel(1 : iSubband);
        isValid = all(candidate >= 0);
        if isValid
            subbandPower(1 : iSubband) = candidate;
        else
            break;
        end
    end
    subbandPower(index) = txPower * (subbandPower ./ sum(subbandPower));

    % * Construct waveform
    infoWaveform = sqrt(transpose(subbandPower)) .* channel' ./ vecnorm(channel, 2, 2)';

    % * Compute capacity
    capacity = sum(log2(1 + subbandPower .* strength / noisePower));

end
