function [fading] = fading_tgn(nSubbands, subbandFrequency, fadingMode)
    % Function:
    %   - simulate channel using the power delay profile of the IEEE TGn NLOS channel model E
    %
    % Input:
    %   - nTxs: number of transmit antennas
    %   - nRxs: number of receive antennas
    %   - nSubbands (N): number of subbands
    %   - subbandFrequency (f_n): the center frequency of subbands
    %   - fadingMode: fading mode "flat" or "selective"
    %
    % Output:
    %   - fading (\boldsymbol{h}) (nTxs * nRxs * nSubbands): fading at each subband
    %
    % Comment:
    %   - for single-user MIMO
    %   - the model only considers power delay profile of clusters
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com) - 07 Mar 20



    % * Define parameters
    nClusters = 4;
    nTaps = 18;
    tapDelay = 1e-9 * [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 560 640 730]';
    tapPower = zeros(nTaps, nClusters);
    tapPower(:, 1) = db2pow([-2.6 -3.0 -3.5 -3.9 -4.5 -5.6 -6.9 -8.2 -9.8 -11.7 -13.9 -16.1 -18.3 -20.5 -22.9 -inf -inf -inf]');
    tapPower(:, 2) = db2pow([-inf -inf -inf -inf -1.8 -3.2 -4.5 -5.8 -7.1 -9.9 -10.3 -14.3 -14.7 -18.7 -19.9 -22.4 -inf -inf]');
    tapPower(:, 3) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -7.9 -9.6 -14.2 -13.8 -18.6 -18.1 -22.8 -inf -inf -inf]');
    tapPower(:, 4) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -20.6 -20.5 -20.7 -24.6]');

    % * Model taps as i.i.d. CSCG variables
    tapGain = repmat(sqrt(tapPower / 2), [1 1 nTxs nUsers]) .* (randn(nTaps, nClusters, nTxs, nUsers) + 1i * randn(nTaps, nClusters, nTxs, nUsers));
    tapGain = permute(sum(tapGain, 2), [1 3 4 2]);
    fading = zeros(nTxs, nSubbands, nUsers);
    switch fadingMode
    case "selective"
        for iTx = 1 : nTxs
            for iSubband = 1 : nSubbands
                for iUser = 1 : nUsers
                    fading(iTx, iSubband, iUser) = sum(tapGain(:, iTx, iUser) .* exp(1i * 2 * pi * subbandFrequency(iSubband) * tapDelay));
                end
            end
        end
    case "flat"
        for iTx = 1 : nTxs
            for iUser = 1 : nUsers
                fading(iTx, :, iUser) = repmat(sum(tapGain(:, iTx, iUser) .* exp(1i * 2 * pi * mean(subbandFrequency) * tapDelay)), [1 nSubbands]);
            end
        end
    end

end
