function [current, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wpt_ff(beta2, beta4, tolerance, directChannel, incidentChannel, reflectiveChannel, irs, txPower, nCandidates, noisePower)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize average output DC current
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - txPower (P): transmit power constraint
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - current (z): harvester output DC current
    %   - irs (\phi) [nReflectors]: IRS reflection coefficients
    %   - infoWavefaaorm (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - initial waveforms by matched filter
    %   - fix splitting ratio and update IRS and waveform iteratively
    %
    % Author & Date: Yang (i@snowztail.com) - 24 Jun 20



    % * Initialize IRS and waveform
    % h^{(0)}, \boldsymbol{M}, \boldsymbol{R}_n
    [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
    % \rho
    powerRatio = 1;
    % \bar{\rho}
    infoRatio = 1 - powerRatio;
    % \bar{R}
    rateConstraint = 0;
    % \boldsymbol{W}_{I/P}^{(0)}
    infoWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
    powerWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
    [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint, tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
    [~, current_] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);

    % * AO
    isConverged = false;
    while ~isConverged
        [irs] = irs_ff(beta2, beta4, nCandidates, rateConstraint, tolerance, infoWaveform, powerWaveform, infoRatio, powerRatio, concatVector, noisePower, concatMatrix, irs);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        [infoWaveform, powerWaveform] = waveform_sdr(beta2, beta4, txPower, nCandidates, rateConstraint, tolerance, infoRatio, powerRatio, noisePower, compositeChannel, infoWaveform, powerWaveform);
        [~, current] = re_sample(beta2, beta4, compositeChannel, noisePower, infoWaveform, powerWaveform, infoRatio, powerRatio);
        isConverged = abs(current - current_) / current <= tolerance || current <= 1e-10;
        current_ = current;
    end

end
