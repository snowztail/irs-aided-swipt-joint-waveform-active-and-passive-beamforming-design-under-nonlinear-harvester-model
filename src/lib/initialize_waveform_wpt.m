function [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wpt(compositeChannel, subbandPower)
    % Function:
    %   - initialize information and power waveform for WPT
    %
    % Input:
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
    %   - optimal WPT, no WIT
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



%     powerRatio = 0.5;
%     infoRatio = 1 - powerRatio;
%     infoWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
%     powerWaveform = sqrt(txPower) * conj(compositeChannel) / norm(compositeChannel);
    powerRatio = 1;
    infoRatio = 1 - powerRatio;
%     infoWaveform = sqrt((1 / 2) * subbandPower) .* exp(1i * angle(conj(compositeChannel)));
    infoWaveform = sqrt(0 * subbandPower) .* exp(1i * angle(conj(compositeChannel)));
    powerWaveform = sqrt(subbandPower) .* exp(1i * angle(conj(compositeChannel)));

end
