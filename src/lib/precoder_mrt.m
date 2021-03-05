function [infoWaveform, powerWaveform] = precoder_mrt(channel, infoAmplitude, powerAmplitude)
    % Function:
    %   - information and power precoder in spatial domain
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude of information waveform in frequency domain
    %   - powerAmplitude (s_P) [1 * nSubbands]: amplitude of power waveform in frequency domain
    %
    % Output:
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %
    % Comment:
    %   - the optimal information and power beamformers coincide at MRT for single-user SWIPT
    %
    % Author & Date: Yang (i@snowztail.com) - 5 Aug 20


	% * Get equivalent channel gain
	channelAmplitude = vecnorm(channel, 2, 2);

    % * Reconstruct waveform by power allocation + beamforming
    infoWaveform = infoAmplitude .* (channel ./ channelAmplitude)';
    powerWaveform = powerAmplitude .* (channel ./ channelAmplitude)';

end
