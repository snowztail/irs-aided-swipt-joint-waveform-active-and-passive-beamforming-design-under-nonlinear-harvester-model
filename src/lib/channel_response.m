function [channel] = channel_response(tapGain, tapDelay, distance, rxGain, subbandFrequency, fadingMode)
    % Function:
    %   - get frequency response of direct, incident and reflective channels
    %
    % Input:
    %   - tapGain [nTaps * nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps * 1]: tap delays
    %   - distance (d): distance between the transmitter and the receiver
    %   - rxGain: receive antenna gain
    %   - subbandFrequency (f_n) [1 * nSubbands]: the center frequency of subbands
    %   - fadingMode: fading mode ('flat', 'selective')
    %
    % Output:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %
    % Comment:
    %   - based on generated tap data
    %   - assume 2 dBi receive antenna gain
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Jun 20


    [pathloss] = path_loss(distance);
    [fading] = multipath_fading(tapGain, tapDelay, subbandFrequency, fadingMode);
    channel = sqrt(pathloss * rxGain) * fading;

end
