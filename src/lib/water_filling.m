function [capacity, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = water_filling(channel, txPower, noisePower)
    % Function:
    %   - water-filling power allocation for OFDM channels
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity: maximum achievable sum rate over all subbands
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
    %
    % Comment:
    %   - for MISO OFDM channels
	%	- convert waveform amplitude to peak value by multiplying sqrt(2)
	%	- assign no power to multisine waveform and set splitting ratio to 0
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Aug 19


    % * Get data
    nSubbands = size(channel, 1);

	% * Obtain SNR and equivalent channel strength
    snr = txPower / noisePower;
	subbandStrength = vecnorm(channel, 2, 2)' .^ 2;

    % * Iterative water-filling power allocation
    subbandPower = zeros(1, nSubbands);
    [subbandStrength, index] = sort(subbandStrength, 'descend');
    baseLevel = 1 ./ (snr * subbandStrength);

    for iSubband = 1 : nSubbands
        waterLevel = (1 + sum(baseLevel(1 : (nSubbands - iSubband + 1)))) / (nSubbands - iSubband + 1);
        subbandPower(1 : (nSubbands - iSubband + 1)) = waterLevel - baseLevel(1 : (nSubbands - iSubband + 1));
		if all(subbandPower >= 0)
			break;
		else
			subbandPower(nSubbands - iSubband + 1) = 0;
		end
    end

    subbandStrength(index) = subbandStrength;
    subbandPower(index) = txPower * (subbandPower ./ sum(subbandPower));

    % * Obtain amplitude and compute capacity
    infoAmplitude = sqrt(2 * subbandPower);
    capacity = sum(log2(1 + subbandPower .* subbandStrength / noisePower));

	% * Initialize multisine waveform and splitting ratio
    powerAmplitude = zeros(1, nSubbands) + eps;
    infoRatio = 1 - eps;
    powerRatio = eps;

end
