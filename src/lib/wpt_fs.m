function [current, infoWaveform, powerWaveform, infoRatio, powerRatio] = wpt_fs(beta2, beta4, tolerance, channel, txPower, nCandidates, noisePower)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize average output DC current
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
    %   - txPower (P): transmit power constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - current (z): harvester output DC current
    %   - infoWavefaaorm (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - no IRS or frequency-selective IRS such that optimal channel is in closed form
    %   - initial waveforms by matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 24 Jun 20



    % * Initialize IRS and waveform
    % \rho
    powerRatio = 1;
    % \bar{\rho}
    infoRatio = 1 - powerRatio;
    % \bar{R}
    rateConstraint = 0;
    % \boldsymbol{W}_{I/P}^{(0)}
    infoWaveform = sqrt(txPower) * conj(channel) / norm(channel);
    powerWaveform = sqrt(txPower) * conj(channel) / norm(channel);
    [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint, tolerance, infoRatio, powerRatio, noisePower, channel, infoWaveform, powerWaveform);
    % z
    [~, current] = re_sample(beta2, beta4, channel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);

end
