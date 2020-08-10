function [tapGain, tapDelay] = tap_tgn(corTx, corRx, losMatrix, variable, propagationMode)
    % Function:
    %   - simulate channel using the power delay profile of the IEEE TGn NLOS channel model E
    %
    % Input:
    %   - corTx (R_t) [nTxs * nTxs]: transmit correlation matrix
    %   - corRx (R_r) [nRxs * nRxs]: receive correlation matrix
    %   - losMatrix [nRxs * nTxs]: fixed LOS response matrix (entries with unit modulus)
    %   - variable {nClusters * nTaps}: i.i.d. CSCG variables as uncorrelated tap response
    %   - propagationMode: propagration mode ('los', 'nlos')
    %
    % Output:
    %   - tapGain [nTaps * nTxs * nRxs]: complex tap gain
    %   - tapDelay [nTaps * 1]: tap delays
    %
    % Comment:
    %   - for single-user MIMO, consider spatial correlation
    %   - only use power delay profile of clusters in the reference (path loss redefined)
    %   - for each tap, the LOS component is fixed while the NLOS component consists of spatially correlated variables
    %   - LOS Ricean factor only apply to the first LOS tap, the remaining taps use NLOS Ricean factor
    %   - LOS matrix depends on geometry position
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com) - 07 Mar 20


    % * Define parameters
    nTxs = size(corTx, 1);
    nRxs = size(corRx, 1);
    nClusters = 4;
    nTaps = 18;
    tapDelay = 1e-9 * [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 560 640 730]';
    tapPower = zeros(nClusters, nTaps);
    tapPower(1, :) = db2pow([-2.6 -3.0 -3.5 -3.9 -4.5 -5.6 -6.9 -8.2 -9.8 -11.7 -13.9 -16.1 -18.3 -20.5 -22.9 -inf -inf -inf]);
    tapPower(2, :) = db2pow([-inf -inf -inf -inf -1.8 -3.2 -4.5 -5.8 -7.1 -9.9 -10.3 -14.3 -14.7 -18.7 -19.9 -22.4 -inf -inf]);
    tapPower(3, :) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -7.9 -9.6 -14.2 -13.8 -18.6 -18.1 -22.8 -inf -inf -inf]);
    tapPower(4, :) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -20.6 -20.5 -20.7 -24.6]);
    losRiceanFactor = db2pow(6);
    nlosRiceanFactor = db2pow(-inf);

    % * Generate tap gains and sum over clusters
    losGain = zeros(nClusters, nTaps, nRxs, nTxs);
    nlosGain = zeros(nClusters, nTaps, nRxs, nTxs);
    tapGain = zeros(nClusters, nTaps, nRxs, nTxs);
    for iCluster = 1 : nClusters
        for iTap = 1 : nTaps
            if iTap == 1 && propagationMode == "los"
                losGain(iCluster, iTap, :, :) = sqrt(losRiceanFactor / (losRiceanFactor + 1)) * losMatrix;
            else
                losGain(iCluster, iTap, :, :) = sqrt(nlosRiceanFactor / (nlosRiceanFactor + 1)) * losMatrix;
            end
            nlosGain(iCluster, iTap, :, :) = sqrt(1 / (nlosRiceanFactor + 1)) * (corRx ^ (1 / 2) * permute(variable(iCluster, iTap, 1 : nRxs, 1 : nTxs), [3 4 1 2]) * corTx ^ (1 / 2));
            tapGain(iCluster, iTap, :, :) = sqrt(tapPower(iCluster, iTap)) * (losGain(iCluster, iTap, :, :) + nlosGain(iCluster, iTap, :, :));
        end
    end

    % * Sum over clusters and reshape to [nTaps, nTxs, nRxs]
    tapGain = permute(sum(tapGain, 1), [2 4 3 1]);

end
