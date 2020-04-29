function [channelAmplitude, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize(txPower, compositeChannel)
    % Function:
    %   - initialize waveform by matched filters
    %
    % InputArg(s):
    %   - txPower: average transmit power
    %   - compositeChannel [\boldsymbol{H}] (nSubbands * 1): total composite channel
    %
    % OutputArg(s):
    %   - channelAmplitude [\boldsymbol{A}] (nSubbands * 1): composite channel amplitude at each subband
    %   - infoAmplitude: [\boldsymbol{A}_I] (nSubbands * 1): amplitude of information weight
    %   - powerAmplitude: [\boldsymbol{A}_P] (nSubbands * 1): amplitude of power weight
    %   - infoRatio: [\bar{\rho}]: information splitting ratio
    %   - powerRatio: \rho: power splitting ratio
    %
    % Comment(s):
    %   - to initialize GP algorithms
    %
    % Author & Date: Yang (i@snowztail.com) - 29 Apr 20



    channelAmplitude = abs(compositeChannel);
    infoAmplitude = sqrt(txPower) * channelAmplitude / norm(channelAmplitude);
    powerAmplitude = sqrt(txPower) * channelAmplitude / norm(channelAmplitude);
    powerRatio = 0.5;
    infoRatio = 1 - powerRatio;

end
