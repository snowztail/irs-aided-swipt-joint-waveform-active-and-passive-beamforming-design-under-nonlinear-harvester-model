function [fading] = fading_tgn(tapGain, tapDelay, nSubbands, subbandFrequency, fadingMode)
    % Function:
    %   - simulate channel using the power delay profile of the IEEE TGn NLOS channel model E
    %
    % Input:
    %   - tapGain [nTaps * nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps]: tap delays
    %   - nSubbands (N): number of subbands
    %   - subbandFrequency (f_n) [nSubbands]: the center frequency of subbands
    %   - fadingMode: fading mode 'flat' or 'selective'
    %
    % Output:
    %   - fading (h) [nSubbands * nTxs * nRxs]: fading at each subband
    %
    % Comment:
    %   - for single-user MIMO
    %   - the model only considers power delay profile of clusters
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com) - 07 Mar 20



    [~, nTxs, nRxs] = size(tapGain);
    fading = zeros(nSubbands, nTxs, nRxs);
    switch fadingMode
    case 'selective'
        for iSubband = 1 : nSubbands
            for iTx = 1 : nTxs
                for iRx = 1 : nRxs
                    fading(iSubband, iTx, iRx) = sum(tapGain(:, iTx, iRx) .* exp(1i * 2 * pi * subbandFrequency(iSubband) * tapDelay));
                end
            end
        end
    case 'flat'
        for iTx = 1 : nTxs
            for iRx = 1 : nRxs
                fading(:, iTx, iRx) = repmat(sum(tapGain(:, iTx, iRx) .* exp(1i * 2 * pi * mean(subbandFrequency) * tapDelay)), [nSubbands 1]);
            end
        end
    end

end
