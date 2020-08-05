function [infoWaveform, powerWaveform] = beamform(channel, infoAmplitude, powerAmplitude)
    % Function:
    %   - obtain waveform in spatial domain
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
    %   - the optimal spatial single-user beamformer is MRT
    %
    % Author & Date: Yang (i@snowztail.com) - 5 Aug 20


    % * Reconstruct waveform by power allocation + beamforming
    infoWaveform = infoAmplitude .* channel' ./ vecnorm(channel, 2, 2)';
    powerWaveform = powerAmplitude .* channel' ./ vecnorm(channel, 2, 2)';

end
