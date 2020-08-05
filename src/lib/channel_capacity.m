function [capacity, infoAmplitude] = channel_capacity(channel, txPower, noisePower)
    % Function:
    %   - calculate the maximum achievable rate based on water-filling power allocation
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity: the maximum achievable rate
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
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
    subbandStrength = transpose(vecnorm(channel, 2, 2)) .^ 2;

    % * Iterative water-filling power allocation (strongest -> weakest)
    subbandPower = zeros(1, nSubbands);
    [subbandStrength, index] = sort(subbandStrength, 'descend');
    baseLevel = 1 ./ (snr * subbandStrength);
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

    % * Restore index and apply power budget
    subbandStrength(index) = subbandStrength;
    subbandPower(index) = txPower * (subbandPower ./ sum(subbandPower));

    % * Obtain amplitude and compute capacity
    infoAmplitude = sqrt(subbandPower);
    capacity = sum(log2(1 + subbandPower .* subbandStrength / noisePower));

end
