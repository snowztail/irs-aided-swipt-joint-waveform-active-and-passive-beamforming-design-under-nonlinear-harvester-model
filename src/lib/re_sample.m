function [rate, current] = re_sample(beta2, beta4, channel, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower)
    % Function:
    %   - compute the output DC current and rate
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - rate (R): user sum rate over all subbands
    %   - current (z): harvester output DC current
    %
    % Comment:
    %   - sample R-E region
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20




    % * Get data
    [nSubbands, nTxs, ~] = size(channel);

    % * Construct coefficient matrices
    % \boldsymbol{W}_{I/P}, \boldsymbol{W}_{I/P,n}
    infoMatrix = vec(infoWaveform) * vec(infoWaveform)';
    powerMatrix = vec(powerWaveform) * vec(powerWaveform)';

    % \boldsymbol{H}_n
    channelMatrix = vec(channel') * vec(channel')';
    channelBlkDiag = block_diagonal(channelMatrix, nTxs, nSubbands);

    % * Compute auxiliaries
    % t'_{I/P,n}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(channelBlkDiag{iSubband + nSubbands} * powerMatrix);
    end
    infoAuxiliary(nSubbands) = hermitianize(infoAuxiliary(nSubbands));
    powerAuxiliary(nSubbands) = hermitianize(powerAuxiliary(nSubbands));

    % * Compute rate and current
    snr = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        snr(iSubband) = infoRatio * abs(channel(iSubband, :) * infoWaveform(:, iSubband)) ^ 2 / noisePower;
    end
    rate = sum(log2(1 + snr));
    current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
        + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
        + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

end
