function [irsPhase, compositeChannel] = irs_fs(directChannel, incidentChannel, reflectiveChannel)
    % Function:
    %   - adjust intelligent reflecting surface to maximize the composite channel strength
    %   - align the extra channel with the direct channel at each subband
    %
    % InputArg(s):
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %
    % OutputArg(s):
    %   - irsPhase (\phi) [nSubbands * nReflectors]: IRS reflection coefficients
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: the composite channel
    %
    % Comment(s):
    %   - assume frequency selective IRS (impractical)
    %   - for single-user scenario
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Apr 20



    [nSubbands, nReflectors, ~] = size(reflectiveChannel);
    compositeChannel = zeros(nSubbands, 1);
    irsPhase = zeros(nSubbands, nReflectors);
    for iSubband = 1 : nSubbands
        for iReflector = 1 : nReflectors
            irsPhase(iSubband, iReflector) = angle(directChannel(iSubband)) - angle(incidentChannel(iSubband, :, iReflector)) - angle(reflectiveChannel(iSubband, iReflector));
        end
        compositeChannel(iSubband) = directChannel(iSubband) + reflectiveChannel(iSubband, :) * diag(exp(1i * irsPhase(iSubband, :))) * squeeze(incidentChannel(iSubband, :, :));
    end

end
