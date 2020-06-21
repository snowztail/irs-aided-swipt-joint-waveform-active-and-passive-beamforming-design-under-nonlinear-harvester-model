function [capacity, infoWaveform, powerWaveform, infoRatio, powerRatio] = wit_no_irs(directChannel, txPower, noisePower)
    % Function:
    %   - optimize the waveform to maximize user rate
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - txPower (P): transmit power constraint
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - capacity (R): channel capacity
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - no IRS
    %   - construct information waveform by water-filling algorithm and matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % p_n^{(0)}
    [capacity, subbandPower] = channel_capacity(directChannel, txPower, noisePower);
    % \rho
    powerRatio = 0;
    % \bar{\rho}
    infoRatio = 1 - powerRatio;
    % \boldsymbol{W}_I^{(0)}
    infoWaveform = sqrt(subbandPower) .* exp(1i * angle(conj(directChannel)));
    % \boldsymbol{W}_P
    powerWaveform = zeros(size(directChannel));

end
