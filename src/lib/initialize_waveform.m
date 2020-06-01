function [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform(txPower, compositeChannel)
    % Function:
    %   - initialize information and power waveform by matched filters
    %
    % Input:
    %   - txPower: average transmit power
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: total composite channel
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - feasible for SISO channel
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    infoWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
    powerWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
    powerRatio = 1;
    infoRatio = 1 - powerRatio;

end
