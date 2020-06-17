function [channel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, distance, tapGain, tapDelay, linkMode)
    % Function:
    %   - get frequency response of direct, incident and reflective channels
    %
    % Input:
    %   - nSubbands (N): number of frequency bands
    %   - subbandFrequency (f_n) [nSubbands]: the center frequency of subbands
    %   - fadingMode: fading mode "flat" or "selective"
    %   - nReflectors: number of reflecting elements in IRS
    %   - distance (d): distance between the transmitter and the receiver
    %   - tapGain [nTaps * nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps]: tap delays
    %   - linkMode: link mode "direct", "incident", or "reflective"
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
    case "incident"
        tapGain = tapGain(:, :, 1 : nReflectors);
    case "reflective"
        tapGain = tapGain(:, 1 : nReflectors, :);
    end

    % * Construct corresponding frequency response
    [pathloss] = path_loss(distance, linkMode);
    [fading] = fading_tgn(tapGain, tapDelay, nSubbands, subbandFrequency, fadingMode);
    channel = fading / sqrt(pathloss);

end
