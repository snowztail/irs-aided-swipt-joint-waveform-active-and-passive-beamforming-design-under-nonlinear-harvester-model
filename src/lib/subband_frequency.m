function [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands)
    % Function:
    %   - calculate the center frequency of subbands
    %
    % Input:
    %   - centerFrequency (f_0): the center frequency between the upper and lower cutoff frequencies
    %   - bandwidth (B): a measure of the width of available frequencies
    %   - nSubbands (N): number of subbands
    %
    % Output:
    %   - subbandFrequency (f_n) [nSubbands]: the center frequency of subbands
    %
    % Comment:
    %   - assume equal spacing between subbands
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20



    subbandFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;

end
