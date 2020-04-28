function [carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands)
    % Function:
    %   - calculate the carrier frequency
    %
    % InputArg(s):
    %   - centerFrequency: a central frequency between the upper and lower cutoff frequencies
    %   - bandwidth: a measure of the width of available frequencies
    %   - nSubbands [N]: number of subbands/subcarriers
    %
    % OutputArg(s):
    %   - carrierFrequency: the frequency of carriers without modulation yet
    %
    % Comment(s):
    %   - assume equal spacing
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20



    carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;

end
