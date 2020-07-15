function [rate, current] = re_sample(beta2, beta4, channel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio)
    % Function:
    %   - compute the output DC current and rate
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - noisePower (\sigma_n^2): average noise power
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Output:
    %   - rate (R): user sum rate over all subbands
    %   - current (z): harvester output DC current
    %
    % Comment:
    %   - sample R-E region
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % * Construct auxiliary variables
    nSubbands = size(infoWaveform, 1);
    % \boldsymbol{H}_{I/P}
    channelMatrix = channel * channel';
    % \boldsymbol{H}_{I/P,n}
    channelCoefMatrix = cell(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        channelCoefMatrix{iSubband + nSubbands} = diag(diag(channelMatrix, iSubband), iSubband);
    end
    % \boldsymbol{W}_{I/P}^{(0)}
    infoMatrix = infoWaveform * infoWaveform';
    powerMatrix = powerWaveform * powerWaveform';
    % t'_{I/P,n}^{(0)}
    infoAuxiliary = zeros(2 * nSubbands - 1, 1);
    powerAuxiliary = zeros(2 * nSubbands - 1, 1);
    for iSubband = - nSubbands + 1 : nSubbands - 1
        infoAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * infoMatrix);
        powerAuxiliary(iSubband + nSubbands) = trace(conj(channelCoefMatrix{iSubband + nSubbands}) * powerMatrix);
    end

    % * Rate
    rate = sum(log2(1 + infoRatio * abs(channel .* infoWaveform) .^ 2 / noisePower));

    % * Current
    current = (1 / 2) * beta2 * powerRatio * (infoAuxiliary(nSubbands) + powerAuxiliary(nSubbands)) ...
        + (3 / 8) * beta4 * powerRatio ^ 2 * (2 * infoAuxiliary(nSubbands) ^ 2 + (powerAuxiliary' * powerAuxiliary)) ...
        + (3 / 2) * beta4 * powerRatio ^ 2 * infoAuxiliary(nSubbands) * powerAuxiliary(nSubbands);

end
