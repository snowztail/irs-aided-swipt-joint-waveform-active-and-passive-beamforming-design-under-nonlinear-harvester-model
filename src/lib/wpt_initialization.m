function [infoWaveform, powerWaveform, infoRatio, powerRatio] = wpt_initialization(channel, txPower)
    % Function:
    %   - initialize the waveform for WPT
    %
    % Input:
    %   - channel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - txPower (P): transmit power constraint
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - no IRS
    %   - construct power waveform by matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Jun 20



    % \rho
    powerRatio = 0.5;
    % \bar{\rho}
    infoRatio = 1 - powerRatio;
    % \boldsymbol{W}_I^{(0)}
    infoWaveform = sqrt(txPower) * conj(channel) / norm(channel);
    % \boldsymbol{W}_P
    powerWaveform = sqrt(txPower) * conj(channel) / norm(channel);

end
