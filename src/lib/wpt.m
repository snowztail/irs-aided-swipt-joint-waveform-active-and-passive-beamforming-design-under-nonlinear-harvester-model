function [current, irs, infoWaveform, powerWaveform, infoRatio, powerRatio] = wpt(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients to maximize average output DC current
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum gain ratio per iteration
    %
    % Output:
    %   - current (z): harvester output DC current
    %   - irs (\phi) [nReflectors * 1]: IRS reflection coefficients
    %   - infoWaveform (w_I) [nTxs * nSubbands]: weight on information waveform
    %   - powerWaveform (w_P) [nTxs * nSubbands]: weight on power waveform
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - solve SDR problem to obtain high-rank IRS outer product matrix
    %   - use Gaussian randomization method to extract IRS vector
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20


    % * Get data
    [nSubbands, nTxs, nReflectors] = size(incidentChannel);

    % * Initialize IRS and composite channel
    % irs = ones(nReflectors, 1);
    irs = exp(1i * 2 * pi * rand(nReflectors, 1));
    [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Initialize waveform and splitting ratio
    infoWaveform = zeros(nTxs, nSubbands) + eps;
    powerWaveform = sqrt(2 * txPower / nSubbands) * compositeChannel' ./ vecnorm(compositeChannel, 2, 2)';
    infoRatio = eps;
    powerRatio = 1 - infoRatio;

    % * AO
    isConverged = false;
    current_ = 0;
    rateConstraint = 0;
    while ~isConverged
        [irs] = irs_sdr(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates, tolerance);
        [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
        % [infoWaveform1, powerWaveform1, ~, ~, ~, current1] = waveform_gp(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio, txPower, noisePower, rateConstraint, tolerance);
        [infoWaveform, powerWaveform, current] = waveform_sdr(beta2, beta4, compositeChannel, infoWaveform, powerWaveform, txPower, nCandidates, tolerance);
        isConverged = abs(current - current_) / current <= tolerance;
        current_ = current;
    end

end
