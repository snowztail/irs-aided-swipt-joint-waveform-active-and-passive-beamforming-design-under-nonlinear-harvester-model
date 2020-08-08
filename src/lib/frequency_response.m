function [channel] = frequency_response(tapGain, tapDelay, distance, nReflectors, subbandFrequency, fadingMode, linkMode)
    % Function:
    %   - get frequency response of direct, incident and reflective channels
    %
    % Input:
    %   - tapGain [nTaps * nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps * 1]: tap delays
    %   - distance (d): distance between the transmitter and the receiver
    %   - nReflectors (L): number of reflecting elements in IRS
    %   - subbandFrequency (f_n) [1 * nSubbands]: the center frequency of subbands
    %   - fadingMode: fading mode 'flat' or 'selective'
    %   - linkMode: link mode 'direct', 'incident', or 'reflective'
    %
    % Output:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %
    % Comment:
    %   - based on generated tap data
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Jun 20


    % * Remove data of unused elements
    switch linkMode
    case 'incident'
        tapGain = tapGain(:, :, 1 : nReflectors);
    case 'reflective'
        tapGain = tapGain(:, 1 : nReflectors, :);
    end

    % * Get pathloss and fading
    [pathloss] = path_loss(distance, linkMode);
    [fading] = fading_tgn(tapGain, tapDelay, subbandFrequency, fadingMode);

    % * Construct frequency response
    channel = sqrt(pathloss) * fading;

end
