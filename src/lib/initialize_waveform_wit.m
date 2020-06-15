function [infoWaveform, powerWaveform, infoRatio, powerRatio] = initialize_waveform_wit(compositeChannel, subbandPower)
    % Function:
    %   - initialize information and power waveform for WIT
    %
    % Input:
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: total composite channel
    %   - subbandPower: optimal power allocation
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %
    % Comment:
    %   - feasible for SISO channel
    %   - optimal WIT, no WPT
    %   - achieve channel capacity
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    powerRatio = 0;
    infoRatio = 1 - powerRatio;
    infoWaveform = sqrt(subbandPower) .* exp(1i * angle(conj(compositeChannel)));
    powerWaveform = zeros(size(compositeChannel));

end
